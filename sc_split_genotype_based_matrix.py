"""
Reference free genotype demultiplexing on pooled scRNA-seq
Jon Xu, Caitlin Falconer
Aug 2018
"""

import pdb
import sys
import numpy as np
import pandas as pd
import math
import datetime
import csv
from scipy.stats import binom


class models:
    def __init__(self, base_calls_mtx, all_POS=[], barcodes=[], num=2, model_genotypes=[], assigned=None, P_c_s=[], P_s_c=[], p_aa_d=[], p_ra_d=[], p_rr_d=[]):
        """
        A complete model class containing SNVs, matrices counts, barcodes, model genotypes with assigned cells 

        Parameters:
             base_calls_mtx: SNV-barcode matrix containing lists of base calls
             all_POS: list of SNVs
             barcodes: list of cell barcodes
             num(int): number of individual model genotypes
             model_genotypes(list(DataFrame): list of num model genotypes represented in snv-probdistrib DataFrame as <RR,RA,AA>
             assigned(list): final lists of cell/barcode assigned to each genotype cluster/model	     
             P_s_c: barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
             P_c_s: barcode/sample matrix containing the probability of seeing barcode c under sample s, whose sum should increase by iterations, as a model indicator
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
        self.model_genotypes = []
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])

        for n in range(self.num):
            self.model_genotypes.append([])
            #self.model_genotypes[n] = pd.DataFrame(np.zeros((len(self.all_POS), 3)),
                #index=all_POS, columns=['RR','RA','AA'])
            # generate random dirichlet distribution to simulate genotypes probability
            self.model_genotypes[n] = pd.DataFrame(np.random.dirichlet((25,50,25),len(self.all_POS)),
                index=self.all_POS, columns=['RR','RA','AA'])

        err = 0.01

        p_d_aa = pd.DataFrame(binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), 1-err),
                index=self.all_POS, columns=self.barcodes)
        p_d_ra = pd.DataFrame(binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), 0.5),
                index=self.all_POS, columns=self.barcodes)
        p_d_rr = pd.DataFrame(binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), err),
                index=self.all_POS, columns=self.barcodes)

        self.p_aa_d = p_d_aa / (p_d_aa + p_d_ra + p_d_rr)
        self.p_ra_d = p_d_ra / (p_d_aa + p_d_ra + p_d_rr)
        self.p_rr_d = p_d_rr / (p_d_aa + p_d_ra + p_d_rr)

    def calculate_model_genotypes(self):
        """
        Update the model genotype based on alternative allele probability for each SNV position based on counts of bases
            Input: self.ref_bc_mtx, self.alt_bc_mtx, P_s_c
            Output: self.model_genotypes
        """

        for n in range(self.num):
            self.model_genotypes[n].loc[:, 'AA'] = self.p_aa_d.dot(self.P_s_c[n]) / (self.p_aa_d.dot(self.P_s_c[n]) + self.p_ra_d.dot(self.P_s_c[n]) + self.p_rr_d.dot(self.P_s_c[n]))
            self.model_genotypes[n].loc[:, 'RA'] = self.p_ra_d.dot(self.P_s_c[n]) / (self.p_aa_d.dot(self.P_s_c[n]) + self.p_ra_d.dot(self.P_s_c[n]) + self.p_rr_d.dot(self.P_s_c[n]))
            self.model_genotypes[n].loc[:, 'RR'] = self.p_rr_d.dot(self.P_s_c[n]) / (self.p_aa_d.dot(self.P_s_c[n]) + self.p_ra_d.dot(self.P_s_c[n]) + self.p_rr_d.dot(self.P_s_c[n]))


    def calculate_cell_likelihood(self):
        """
        Calculate cell|sample likelihood and derive sample|cell probability
        Input: self.model_genotypes, self.ref_bc_mtx, self.alt_bc_mtx
        Output: self.P_c_s, self.P_s_c
        """
        err = 0.01

        for n in range(self.num):
            rr = pd.DataFrame(data=binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), err), index=self.alt_bc_mtx.index, columns=self.alt_bc_mtx.columns)
            ra = pd.DataFrame(data=binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), 0.5), index=self.alt_bc_mtx.index, columns=self.alt_bc_mtx.columns)
            aa = pd.DataFrame(data=binom.pmf(self.alt_bc_mtx, (self.alt_bc_mtx + self.ref_bc_mtx), 1-err), index=self.alt_bc_mtx.index, columns=self.alt_bc_mtx.columns)
            matcalc = rr.multiply(self.model_genotypes[n].loc[:,'RR'], axis=0) + \
                      ra.multiply(self.model_genotypes[n].loc[:,'RA'], axis=0) + \
                      aa.multiply(self.model_genotypes[n].loc[:,'AA'], axis=0)
            
            #matcalc = pd.DataFrame(np.zeros_like(self.ref_bc_mtx), index=self.ref_bc_mtx.index, columns=self.ref_bc_mtx.columns)
            self.P_c_s.loc[:, n] = 2 ** (matcalc.apply(np.log2).sum(axis=0))
        self.P_s_c = self.P_c_s.div(self.P_c_s.sum(axis=1), axis=0)


    def assign_cells(self):
        """
	    final assignment of cells according to model genotypes
	    """

        for barcode in self.barcodes:  # for each of all cells

            # if likelihood is equal in cases such as: barcode has no coverage over any snv region
            # actual cell assignment step
            if self.P_s_c.loc[barcode,:].max() >= 0.8:
                n = pd.Series.idxmax(self.P_s_c.loc[barcode,:])
                self.assigned[n].append(barcode)

        for n in range(self.num):
            self.assigned[n] = sorted(self.assigned[n])


def run_model(base_calls_mtx, num_models):

    model = models(base_calls_mtx, num_models)

    print("Commencing E-M")
    
    iterations = 0
    sum_log_likelihoods = []

    while iterations < 6:

        iterations += 1

        print("calculating cell likelihood ", datetime.datetime.now().time())
        model.calculate_cell_likelihood()

        print("calculating model ", datetime.datetime.now().time())
        model.calculate_model_genotypes()

        sum_log_likelihood = model.P_c_s.apply(np.log10).sum().sum()
        sum_log_likelihoods.append(sum_log_likelihood)
        print("log likelihood of iteration {}".format(iterations), sum_log_likelihood)

    model.assign_cells()

    for n in range(num_models):
        with open('barcodes_{}_f_simple_6.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')

        model.model_genotypes[n].to_csv('model_genotypes{}_f_simple_6.csv'.format(n))

    model.P_c_s.to_csv('P_c_s_f_simple_6.csv')
    model.P_s_c.to_csv('P_s_c_f_simple_6.csv')
    
    print("Finished model at {}".format(datetime.datetime.now().time()))
    print(sum_log_likelihoods)


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
    # ref_csv = 'ref_filtered.csv'  # reference matrix
    # alt_csv = 'alt_filtered.csv'  # alternative matrix

    ref_csv = 'test_ref.csv'
    alt_csv = 'test_alt.csv'

    print ("Starting data collection", datetime.datetime.now().time())
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, num_models)

if __name__ == "__main__":
    main()

