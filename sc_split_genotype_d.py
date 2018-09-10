"""
Reference free genotype demultiplexing on pooled scRNA-seq
Jon Xu
Lachlan Coin
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
             lP_c_s(DataFrame: barcode/sample matrix containing the log likelihood of seeing barcode c under sample s, whose sum should increase by each iteration
             assigned(list): final lists of cell/barcode assigned to each genotype cluster/model
             model_genotypes(list(DataFrame): list of num model genotypes represented in snv-prob distrib DataFrame as <RR,RA,AA>
             p_d_aa, p_d_ra, p_d_rr(DataFrame): SNV/barcode matrix containing P(D|AA), P(D|RA), P(D|RR)
             p_aa_d, p_ra_d, p_rr_d(DataFrame): SNV/barcode matrix containing P(AA|D), P(RA|D), P(RR|D)

        """

        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.all_POS = self.ref_bc_mtx.index.values.tolist()
        self.barcodes = self.ref_bc_mtx.columns.values.tolist()
        self.num = num + 1  # including an additional background state for doublets
        self.P_s_c = pd.DataFrame(np.zeros((len(self.barcodes), self.num)),
                    index = self.barcodes, columns = list(range(self.num)))
        self.lP_c_s = pd.DataFrame(np.zeros((len(self.barcodes), self.num)),
                    index = self.barcodes, columns = list(range(self.num)))
        self.model_genotypes = []
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])

        self.model_genotypes.append([])
        # set background genotype for each SNVs in the model
        self.model_genotypes[0] = pd.DataFrame(np.tile([0.25,0.5,0.25],[len(self.all_POS),1]),
                index=self.all_POS, columns=['RR','RA','AA'])
        for n in range(1, self.num):
            self.model_genotypes.append([])
            # generate random dirichlet distribution to simulate genotypes probability
            self.model_genotypes[n] = pd.DataFrame(np.random.dirichlet((25,50,25),len(self.all_POS)),
                index=self.all_POS, columns=['RR','RA','AA'])

        err = 0.01  # error rate assumption

        # binomial probability for P(D|AA,RA,RR) with the alt count vs total count condition and (err, 0.5, 1-err) genotype probability
        self.p_d_aa = pd.DataFrame(binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), 1-err), index=self.all_POS, columns=self.barcodes)
        self.p_d_ra = pd.DataFrame(binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), 0.5), index=self.all_POS, columns=self.barcodes)
        self.p_d_rr = pd.DataFrame(binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), err), index=self.all_POS, columns=self.barcodes)

        # transform into posterior probability P(AA,RA,RR|D)
        self.p_aa_d = self.p_d_aa / (self.p_d_aa + self.p_d_ra + self.p_d_rr)
        self.p_ra_d = self.p_d_ra / (self.p_d_aa + self.p_d_ra + self.p_d_rr)
        self.p_rr_d = self.p_d_rr / (self.p_d_aa + self.p_d_ra + self.p_d_rr)


    def calculate_model_genotypes(self):
        """
        Update the model genotype by simulating the count distribution using P(s|c) and P(g|D) of each barcode on a certain snv to the model

        """

        for n in range(self.num):
            self.model_genotypes[n].loc[:, 'AA'] = self.p_aa_d.dot(self.P_s_c[n]) / (self.p_aa_d.dot(self.P_s_c[n]) + self.p_ra_d.dot(self.P_s_c[n]) + self.p_rr_d.dot(self.P_s_c[n]))
            self.model_genotypes[n].loc[:, 'RA'] = self.p_ra_d.dot(self.P_s_c[n]) / (self.p_aa_d.dot(self.P_s_c[n]) + self.p_ra_d.dot(self.P_s_c[n]) + self.p_rr_d.dot(self.P_s_c[n]))
            self.model_genotypes[n].loc[:, 'RR'] = self.p_rr_d.dot(self.P_s_c[n]) / (self.p_aa_d.dot(self.P_s_c[n]) + self.p_ra_d.dot(self.P_s_c[n]) + self.p_rr_d.dot(self.P_s_c[n]))
        self.model_genotypes[0] = pd.DataFrame(np.tile([0.25,0.5,0.25],[len(self.all_POS),1]),index=self.all_POS, columns=['RR','RA','AA'])   # reset the background genotype

    def calculate_cell_likelihood(self):
        """
        Calculate cell|sample likelihood P(c|s) and derive sample|cell probability P(s|c)
        P(c|s_v) = P(D|AA) * P(g_AA|s) + P(D|RA) * P(g_RA|s) + P(D|RR) * P(g_RR|s)
        log(P(c|s)) = sum_v(log(P(c|s_v)))
        P(s_n|c) = P(c|s_n) / [P(c|s_1) + P(c|s_2) + ... + P(c|s_n)]

        """

        P_s = []

        for n in range(self.num):
            matcalc = self.p_d_rr.multiply(self.model_genotypes[n].loc[:,'RR'], axis=0) + \
                      self.p_d_ra.multiply(self.model_genotypes[n].loc[:,'RA'], axis=0) + \
                      self.p_d_aa.multiply(self.model_genotypes[n].loc[:,'AA'], axis=0)            
            self.lP_c_s.loc[:, n] = matcalc.apply(np.log2).sum(axis=0)  # log likelihood to avoid python computation limit of 2^+/-308
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
    while iterations < 15:

        iterations += 1
        print("Iteration {}".format(iterations))
        model.calculate_cell_likelihood()
        print("cell origin probabilities ", model.P_s_c)
        model.calculate_model_genotypes()
        print("model_genotype_d: ", model.model_genotypes)
        sum_log_likelihood.append(model.lP_c_s.sum().sum())

    model.assign_cells()

    # generate outputs
    for n in range(num_models+1):
        with open('barcodes_genotype_d_{}.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')    
    model.P_s_c.to_csv('P_s_c_genotype_d.csv')
    print(sum_log_likelihood)
    print("Finished model at {}".format(datetime.datetime.now().time()))


def read_base_calls_matrix(ref_csv, alt_csv):

    """ Read in an existing matrix from a csv file"""
    
    base_calls_mtx = []
    print('reading in reference matrix')
    base_calls_mtx.append(pd.read_csv(ref_csv, header=0, index_col=0))
    print('reading in alternate matrix')
    base_calls_mtx.append(pd.read_csv(alt_csv, header=0, index_col=0))
    print("Base call matrix finished", datetime.datetime.now().time())
    return base_calls_mtx


def main():

    num_models = 2          # number Fof models in each run

    # Mixed donor files
    ref_csv = 'ref_filtered.csv'  # reference matrix
    alt_csv = 'alt_filtered.csv'  # alternative matrix

    print ("Starting data collection", datetime.datetime.now().time())
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, num_models)

if __name__ == "__main__":
    main()

