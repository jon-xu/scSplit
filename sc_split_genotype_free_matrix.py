"""
Reference free genotype demultiplexing on pooled scRNA-seq
Jon Xu, Caitlin Falconer
Aug 2018
"""

import pdb
import sys
import math
import numpy as np
import pandas as pd
from scipy.stats import binom
import datetime
import csv

class models:
    def __init__(self, base_calls_mtx, num=2):
        """
        A complete model class containing SNVs, matrices counts, barcodes, model genotypes with assigned cells 

        Parameters:
             base_calls_mtx(list of Dataframes): SNV-barcode matrix containing lists of base calls
             all_POS(list): list of SNVs positions
             barcodes(list): list of cell barcodes
             num(int): number of total samples
             P_s_c(DataFrame): barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
             P_c_s(DataFrame: barcode/sample matrix containing the likelihood of seeing barcode c under sample s, whose sum should increase by each iteration
             assigned(list): final lists of cell/barcode assigned to each genotype cluster/model
             model_genotypes(list): list of num model genotypes based on P(AA)

        """

        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.all_POS = self.ref_bc_mtx.index.values.tolist()
        self.barcodes = self.ref_bc_mtx.columns.values.tolist()
        self.num = num
        self.P_s_c = pd.DataFrame(np.zeros((len(self.barcodes), self.num)),
                    index = self.barcodes, columns = list(range(self.num)))
        self.P_c_s = pd.DataFrame(np.zeros((len(self.barcodes), self.num)),
                    index = self.barcodes, columns = list(range(self.num)))
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])
        self.model_genotypes = pd.DataFrame(np.zeros((len(self.all_POS), self.num)),
                    index=self.all_POS, columns=range(self.num))
        for n in range(self.num):
            beta_sim = np.random.beta(self.ref_bc_mtx.sum().sum(), self.alt_bc_mtx.sum().sum(), size = (len(self.all_POS), self.num))
            self.model_genotypes.loc[:, n] = [item[1] for item in beta_sim]

    def calculate_model_genotypes(self):
        """
        Update the model genotype by distributing the alt and total counts of each barcode on a certain snv to the model based on P(s|c)

        """

        self.model_genotypes = (self.alt_bc_mtx.dot(self.P_s_c) + 1) / ((self.alt_bc_mtx + self.ref_bc_mtx).dot(self.P_s_c) + 2)


    def calculate_cell_likelihood(self):
        """
        Calculate cell|sample likelihood P(c|s) and derive sample|cell probability P(s|c)
        P(c|s_v) = P(N(A),N(R)|s) = P(g_A|s)^N(A) * (1-P(g_A|s))^N(R)
        log(P(c|s)) = sum_v{(N(A)_c,v*log(P(g_A|s)) + N(R)_c,v*log(1-P(g_A|s)))}
        P(s_n|c) = P(c|s_n) / [P(c|s_1) + P(c|s_2) + ... + P(c|s_n)]

        """

        for n in range(self.num):
            matcalc = self.alt_bc_mtx.multiply(self.model_genotypes.loc[:, n].apply(np.log2), axis=0) \
                    + self.ref_bc_mtx.multiply((1 - self.model_genotypes.loc[:, n]).apply(np.log2), axis=0)
            self.P_c_s.loc[:, n] = 2 ** (matcalc.sum(axis=0) + 800)
        self.P_s_c = self.P_c_s.div(self.P_c_s.sum(axis=1), axis=0)


    def assign_cells(self):
        """
	    Final assignment of cells according to P(s|c) > 0.8

	    """

        for n in range(self.num):
            self.assigned[n] = sorted(self.P_s_c.loc[self.P_s_c[n] >= 0.8])


def run_model(base_calls_mtx, num_models):

    model = models(base_calls_mtx, num_models)

    print("Commencing E-M")
    
    iterations = 0
    sum_log_likelihoods = []

    while iterations < 3:
        iterations += 1
        print("Iteration {}".format(iterations))

        print("calculating cell likelihood ", datetime.datetime.now().time())
        model.calculate_cell_likelihood()
        print("cell origin probabilities ", model.P_s_c)

        print("calculating model ", datetime.datetime.now().time())
        model.calculate_model_genotypes()
        print("model_genotypes: ", model.model_genotypes)

        sum_log_likelihood = model.P_c_s.apply(np.log10).sum().sum()
        sum_log_likelihoods.append(sum_log_likelihood)
        print("log likelihood of iteration {}".format(iterations), sum_log_likelihood)

    model.assign_cells()

    for n in range(num_models):
        with open('barcodes_{}_f_simple_3.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')
    
    print("Finished model at {}".format(datetime.datetime.now().time()))
    print(sum_log_likelihoods)


def read_base_calls_matrix(ref_csv, alt_csv):

    """ Read in existing matrix from the csv files """

    base_calls_mtx = []
    print('reading in reference matrix')
    base_calls_mtx.append(pd.read_csv(ref_csv, header=0, index_col=0))
    print('reading in alternate matrix')
    base_calls_mtx.append(pd.read_csv(alt_csv, header=0, index_col=0))
    print("Base call matrix finished", datetime.datetime.now().time())
    return base_calls_mtx


def main():

    num_models = 2          # number of models in each run

    # Mixed donor files
    # bc_file = "bc_sorted.txt"   # validated cell barcodes
    # ref_csv = 'ref_filtered.csv'  # reference matrix
    # alt_csv = 'alt_filtered.csv'  # alternative matrix

    bc_file = 'test.txt'
    ref_csv = 'test_ref.csv'
    alt_csv = 'test_alt.csv'

    print ("Starting data collection", datetime.datetime.now().time())
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, num_models)

if __name__ == "__main__":
    main()

