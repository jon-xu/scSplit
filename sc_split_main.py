"""
Reference free MAF-based demultiplexing on pooled scRNA-seq
Jon Xu (jun.xu@uq.edu.au)
Lachlan Coin
Aug 2018
"""

import sys
import math
import numpy as np
import pandas as pd
from scipy.stats import binom
from scipy.sparse import csr_matrix
import datetime
import csv

class models:
    def __init__(self, base_calls_mtx, num=2):
        """
        A complete model class containing SNVs, matrices counts, barcodes, model minor allele frequency with assigned cells 

        Parameters:
             base_calls_mtx(list of Dataframes): SNV-barcode matrix containing lists of base calls
             all_POS(list): list of SNVs positions
             barcodes(list): list of cell barcodes
             num(int): number of total samples
             P_s_c(DataFrame): barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
             lP_c_s(DataFrame): barcode/sample matrix containing the log likelihood of seeing barcode c under sample s, whose sum should increase by each iteration
             assigned(list): final lists of cell/barcode assigned to each cluster/model
             model_MAF(list): list of num model minor allele frequencies based on P(A)

        """

        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.all_POS = base_calls_mtx[2].tolist()
        self.barcodes = base_calls_mtx[3].tolist()
        self.num = num + 1  # including an additional background state for doublets
        self.P_s_c = pd.DataFrame(np.zeros((len(self.barcodes), self.num)), index = self.barcodes, columns = range(self.num))
        self.lP_c_s = pd.DataFrame(np.zeros((len(self.barcodes), self.num)), index = self.barcodes, columns = range(self.num))
        self.A_s_c = pd.DataFrame(np.zeros((len(self.barcodes), self.num)), index = self.barcodes, columns = range(self.num))
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])
        self.model_MAF = pd.DataFrame(np.zeros((len(self.all_POS), self.num)), index=self.all_POS, columns=range(self.num))
        # set background alt count proportion as fixed minor allele frequency for each SNVs in the model, pseudo count is added for 0 counts on multi-base SNPs
        self.model_MAF.loc[:, 0] = (self.alt_bc_mtx.sum(axis=1) + 1) / (self.ref_bc_mtx.sum(axis=1) + self.alt_bc_mtx.sum(axis=1) + 2)
        for n in range(1, self.num):
            # use total ref count and alt count to generate probability simulation
            beta_sim = np.random.beta(self.ref_bc_mtx.sum(), self.alt_bc_mtx.sum(), size = (len(self.all_POS), 1))
            self.model_MAF.loc[:, n] = [1 - item[0] for item in beta_sim]   # P(A) = 1 - P(R)


    def calculate_model_MAF(self):
        """
        Update the model minor allele frequency by distributing the alt and total counts of each barcode on a certain snv to the model based on P(s|c)

        """

        self.model_MAF = pd.DataFrame((self.alt_bc_mtx.dot(self.P_s_c) + 1) / ((self.alt_bc_mtx + self.ref_bc_mtx).dot(self.P_s_c) + 2),
                                        index = self.all_POS, columns = range(self.num))
        self.model_MAF.loc[:, 0] = self.model_MAF.loc[:, 1:(self.num-1)].mean(axis=1)   # reset the background MAF


    def calculate_cell_likelihood(self):
        """
        Calculate cell|sample likelihood P(c|s) and derive sample|cell probability P(s|c)
        P(c|s_v) = P(N(A),N(R)|s) = P(g_A|s)^N(A) * (1-P(g_A|s))^N(R)
        log(P(c|s)) = sum_v{(N(A)_c,v*log(P(g_A|s)) + N(R)_c,v*log(1-P(g_A|s)))}
        P(s_n|c) = P(c|s_n) / [P(c|s_1) + P(c|s_2) + ... + P(c|s_n)]

        """

        P_s = []

        for n in range(self.num):
            matcalc = self.alt_bc_mtx.T.multiply(self.model_MAF.loc[:, n].apply(np.log2)).T \
                    + self.ref_bc_mtx.T.multiply((1 - self.model_MAF.loc[:, n]).apply(np.log2)).T
            self.lP_c_s.loc[:, n] = matcalc.sum(axis=0).tolist()[0]  # log likelihood to avoid python computation limit of 1e-323/1e+308
            if n == 0:
                P_s.append(0.02)  # probability of doublet ratio
            else:
                P_s.append((1 - 0.02) / (self.num - 1))  # even distribution of P(s) across all other singlet samples
    
        # log(P(s1|c) = log{1/[1+P(c|s2)/P(c|s1)]} = -log[1+P(c|s2)/P(c|s1)] = -log[1+2^(logP(c|s2)-logP(c|s1))]
        for i in range(self.num):
            denom = 0
            for j in range(self.num):
                denom += 2 ** (self.lP_c_s.loc[:, j] + np.log2(P_s[j]) - self.lP_c_s.loc[:, i] - np.log2(P_s[i]))
            self.P_s_c.loc[:, i] = 1 / denom


    def assign_cells(self):
        """
	    Final assignment of cells according to P(s|c) >= 0.9

	    """

        for n in range(self.num):
            self.assigned[n] = sorted(self.P_s_c.loc[self.P_s_c[n] >= 0.9].index.values.tolist())



def run_model(base_calls_mtx, num_models):

    model = models(base_calls_mtx, num_models)
    
    iterations = 0
    sum_log_likelihood = []

    # commencing E-M
    while iterations < 10:
        iterations += 1
        progress = 'Iteration ' + str(iterations) + '   ' + str(datetime.datetime.now()) + '\n'
        with open('wip.log', 'a') as myfile:
            myfile.write(progress)
        model.calculate_cell_likelihood()  # E-step, calculate the expected cell origin likelihood with a function of model.model_MAF (theta)
        model.calculate_model_MAF()  # M-step, to optimise unknown model parameter model.model_MAF (theta)
        # approximation due to python calculation limit
        sum_log_likelihood.append(model.lP_c_s.max(axis=1).sum())  # L = Sum_c{log(Sum_s(P(c|s))}
        # sum_log_likelihood.append(((2**model.lP_c_s).sum(axis=1)+1e-323).apply(np.log2).sum())  # L = Sum_c{log(Sum_s(P(c|s))}

    model.assign_cells()

    # generate outputs
    for n in range(num_models+1):
        with open('barcodes_{}.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')
    model.P_s_c.to_csv('P_s_c.csv')
    model.model_MAF.to_csv('model_maf.csv')
    print(sum_log_likelihood)
    print("Finished model at {}".format(datetime.datetime.now().time()))


def read_base_calls_matrix(ref_csv, alt_csv):

    """ Read in existing matrix from the csv files """

    ref = pd.read_csv(ref_csv, header=0, index_col=0)  # read ref matrix with header line and column
    alt = pd.read_csv(alt_csv, header=0, index_col=0)  # read ref matrix with header line and column
    ref_s = csr_matrix(ref.values)
    alt_s = csr_matrix(alt.values)
    base_calls_mtx = [ref_s, alt_s, ref.index, ref.columns]
    print("Base call matrix finished", datetime.datetime.now().time())
    return base_calls_mtx


def main():

    num_models = 2          # number of models in each run

    # input and output files
    ref_csv = 'ref_filtered.csv'  # reference matrix
    alt_csv = 'alt_filtered.csv'  # alternative matrix

    print ("Starting data collection", datetime.datetime.now().time())
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, num_models)

if __name__ == "__main__":
    main()

