"""
Reference free AF-based demultiplexing on pooled scRNA-seq
Jon Xu (jun.xu@uq.edu.au)
Lachlan Coin
Aug 2018
"""

import sys
import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import datetime
import csv

class models:
    def __init__(self, base_calls_mtx, seed_vcf, num):
        """
        A complete model class containing SNVs, matrices counts, barcodes, model allele fraction with assigned cells 

        Parameters:
             base_calls_mtx(list of Dataframes): SNV-barcode matrix containing lists of base calls
             all_POS(list): list of SNVs positions
             barcodes(list): list of cell barcodes
             num(int): number of total samples
             P_s_c(DataFrame): barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
             lP_c_s(DataFrame): barcode/sample matrix containing the log likelihood of seeing barcode c under sample s, whose sum should increase by each iteration
             assigned(list): final lists of cell/barcode assigned to each cluster/model
             model_af(list): list of num model allele frequencies based on P(A)

        """

        err = 0.01  # error rate for model seeds
        dbl = 0.02  # doublet ratio assumption
        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.all_POS = base_calls_mtx[2].tolist()
        self.barcodes = base_calls_mtx[3].tolist()
        self.num = num + 1  # including an additional background state for doublets
        self.P_s_c = pd.DataFrame(0, index = self.barcodes, columns = range(self.num))
        self.lP_c_s = pd.DataFrame(0, index = self.barcodes, columns = range(self.num))
        self.assigned = []
        self.P_s = [dbl]  # assuming P_s[0], i.e doublet has 2% probability
        for _ in range(self.num):
            self.assigned.append([])
        self.model_af = pd.DataFrame(0, index=self.all_POS, columns=range(self.num))
        # set background alt count proportion as fixed allele fraction for each SNVs in the model, pseudo count is added for 0 counts on multi-base SNPs
        self.model_af.loc[:, 0] = (self.alt_bc_mtx.sum(axis=1) + 1) / (self.ref_bc_mtx.sum(axis=1) + self.alt_bc_mtx.sum(axis=1) + 2)
        for n in range(1, self.num):
            self.P_s.append((1 - dbl) / (self.num - 1))  # even initial distribution of P(s) across all other singlet samples
            # use total ref count and alt count on each position of csr_matrix to generate probability simulation using beta distribution
            N_A = self.alt_bc_mtx.sum(axis=1) + 1
            N_R = self.ref_bc_mtx.sum(axis=1) + 1
            N_T = N_A + N_R
            self.model_af.loc[:, n] = [item[0] for item in np.random.beta(100*N_A/N_T, 100*N_R/N_T)]

        last = ''
        in_vcf = vcf.Reader(open(seed_vcf, 'r'))
        for record in in_vcf:
            if (record.ID != last) & (record.ID in self.all_POS):   # only need those SNVs captured in our matrices
                for n in range(1, self.num):
                    self.model_af.loc[record.ID, n] = 0.5 * record.samples[n-1]['GP'][1] + record.samples[n-1]['GP'][2]
            last = record.ID

        self.model_af[self.model_af == 0] = err
        self.model_af[self.model_af == 1] = 1 - err


    def calculate_model_af(self):
        """
        Update the model allele fraction by distributing the alt and total counts of each barcode on a certain snv to the model based on P(s|c)

        """

        pseudo_count = 0.01  # pseudo count for zero entries
        N_ref = self.ref_bc_mtx.sum(axis=1)
        N_alt = self.alt_bc_mtx.sum(axis=1)
        N_ref[N_ref == 0] = pseudo_count
        N_alt[N_alt == 0] = pseudo_count
        k_ref = N_ref / (N_ref + N_alt)
        k_alt = N_alt / (N_ref + N_alt)
        self.model_af = pd.DataFrame((self.alt_bc_mtx.dot(self.P_s_c) + k_alt) / ((self.alt_bc_mtx + self.ref_bc_mtx).dot(self.P_s_c) + k_ref + k_alt),
                                        index = self.all_POS, columns = range(self.num))
        self.model_af.loc[:, 0] = self.model_af.loc[:, 1:(self.num-1)].mean(axis=1)   # reset the background AF


    def calculate_cell_likelihood(self):
        """
        Calculate cell|sample likelihood P(c|s) and derive sample|cell probability P(s|c)
        P(c|s_v) = P(N(A),N(R)|s) = P(g_A|s)^N(A) * (1-P(g_A|s))^N(R)
        log(P(c|s)) = sum_v{(N(A)_c,v*log(P(g_A|s)) + N(R)_c,v*log(1-P(g_A|s)))}
        P(s_n|c) = P(c|s_n) / [P(c|s_1) + P(c|s_2) + ... + P(c|s_n)]

        """

        for n in range(self.num):
            matcalc = self.alt_bc_mtx.T.multiply(self.model_af.loc[:, n].apply(np.log2)).T \
                    + self.ref_bc_mtx.T.multiply((1 - self.model_af.loc[:, n]).apply(np.log2)).T
            self.lP_c_s.loc[:, n] = matcalc.sum(axis=0).tolist()[0]  # log likelihood to avoid python computation limit of 1e-323/1e+308
    
        # log(P(s1|c) = log{1/[1+P(c|s2)/P(c|s1)]} = -log[1+P(c|s2)/P(c|s1)] = -log[1+2^(logP(c|s2)-logP(c|s1))]
        for i in range(self.num):
            denom = 0
            for j in range(self.num):
                denom += 2 ** (self.lP_c_s.loc[:, j] + np.log2(self.P_s[j]) - self.lP_c_s.loc[:, i] - np.log2(self.P_s[i]))
            self.P_s_c.loc[:, i] = 1 / denom


    def assign_cells(self):
        """
	    Final assignment of cells according to P(s|c) >= 0.9

	    """

        for n in range(self.num):
            self.assigned[n] = sorted(self.P_s_c.loc[self.P_s_c[n] >= 0.9].index.values.tolist())



def run_model(base_calls_mtx, seed_vcf, num_models):

    model = models(base_calls_mtx, seed_vcf, num_models)
    
    iterations = 0
    sum_log_likelihood = [1,2]  # dummy likelihood as a start

    # commencing E-M
    while sum_log_likelihood[-2] != sum_log_likelihood[-1]:
        iterations += 1
        progress = 'Iteration ' + str(iterations) + '   ' + str(datetime.datetime.now()) + '\n'
        with open('wip.log', 'a') as myfile: myfile.write(progress)
        model.calculate_cell_likelihood()  # E-step, calculate the expected cell origin likelihood with a function of model.model_af (theta)
        model.calculate_model_af()  # M-step, to optimise unknown model parameter model.model_af (theta)
        model.model_af.to_csv('model_af_' + str(iterations))
        # approximation due to python calculation limit
        sum_log_likelihood.append(model.lP_c_s.max(axis=1).sum())  # L = Prod_c[Sum_s(P(c|s))], thus LL = Sum_c{log[Sum_s(P(c|s))]}
        # sum_log_likelihood.append(((2**model.lP_c_s).sum(axis=1)+1e-323).apply(np.log2).sum())

    model.assign_cells()

    # generate outputs
    for n in range(num_models+1):
        with open('barcodes_{}.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')
    model.P_s_c.to_csv('P_s_c.csv')
    model.model_af.to_csv('model_af.csv')
    print(sum_log_likelihood)
    progress = 'scSplit finished at: ' + str(datetime.datetime.now()) + '\n'
    with open('wip.log', 'a') as myfile: myfile.write(progress)


def read_base_calls_matrix(ref_csv, alt_csv):

    """ Read in existing matrix from the csv files """

    ref = pd.read_csv(ref_csv, header=0, index_col=0)  # read ref matrix with header line and column
    alt = pd.read_csv(alt_csv, header=0, index_col=0)  # read alt matrix with header line and column
    ref_s = csr_matrix(ref.values)
    alt_s = csr_matrix(alt.values)
    base_calls_mtx = [ref_s, alt_s, ref.index, ref.columns]
    progress = 'AF matrices uploaded: ' + str(datetime.datetime.now()) + '\n'
    with open('wip.log', 'a') as myfile: myfile.write(progress)
    return base_calls_mtx


def main():

    num_models = 8          # number of models in each run

    # input and output files
    ref_csv = 'ref_filtered.csv'  # reference matrix
    alt_csv = 'alt_filtered.csv'  # alternative matrix
    seed_vcf = 'input8.vcf'  # vcf serving as model seed

    progress = 'Starting data collection: ' + str(datetime.datetime.now()) + '\n'
    with open('wip.log', 'a') as myfile: myfile.write(progress)
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, seed_vcf, num_models)

if __name__ == '__main__':
    main()


