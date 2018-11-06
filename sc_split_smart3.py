"""
Reference free AF-based demultiplexing on pooled scRNA-seq (second intialisation, but including doublets and excluding multiplet)
Jon Xu (jun.xu@uq.edu.au)
Lachlan Coin
Aug 2018
"""

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import datetime
import csv

class models:
    def __init__(self, base_calls_mtx, num=2):
        """
        Model class containing SNVs, matrices counts, barcodes, model allele fraction with assigned cells 

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
        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.all_POS = base_calls_mtx[2].tolist()
        self.barcodes = base_calls_mtx[3].tolist()
        self.num = num + int(num * (num -1) / 2) + 1  # including additional background states for doublets and one for all other multiplets
        self.singlets = num  # number of singlet states
        self.P_s_c = pd.DataFrame(0, index = self.barcodes, columns = range(self.num))
        self.lP_c_s = pd.DataFrame(0, index = self.barcodes, columns = range(self.num))
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])
        self.model_af = pd.DataFrame(0, index=self.all_POS, columns=range(self.num))
        self.seeds = [np.argmax((self.ref_bc_mtx + self.alt_bc_mtx).sum(axis=0))]
        self.pseudo = 0.01


    def initialise_model(self, x):
        
        ### initialise multiplet:

        dbl = 0.02  # doublet ratio assumption
        # allele counts with pseudo count on each SNV position
        N_A = self.alt_bc_mtx.sum(axis=1) + self.pseudo
        N_R = self.ref_bc_mtx.sum(axis=1) + self.pseudo
        N_T = N_A + N_R
        k_ref = N_R / N_T
        k_alt = N_A / N_T
        # set background alt count proportion as allele fraction for each SNVs of doublet state, with pseudo count added for 0 counts on multi-base SNPs
        self.model_af.loc[:, 0] = N_A / N_T
        self.P_s = []  # assuming P_s[0], i.e doublet has 2% probability

        ### initialise singlets:

        if x <= self.singlets:  # initialise with selected seeds
            # initialise all states using two different strategies (cell info for n <= x, all cell counts for n > x)
            for n in range(1, self.singlets + 1):    
                self.P_s.append((1 - dbl) / (self.singlets))  # even initial distribution of P(s) across all other singlet samples
                if n <= x:  # use seeded cells to initialise seeded states the model
                    barcode_alt = self.alt_bc_mtx.getcol(self.seeds[n-1]).toarray()
                    barcode_ref = self.ref_bc_mtx.getcol(self.seeds[n-1]).toarray()
                    self.model_af.loc[:, n] = (barcode_alt + k_alt) / (barcode_alt + barcode_ref + k_alt + k_ref)
                else:
                    # for the rest states, use total ref/alt count on each SNV position within sparse matrices to generate probability simulation using beta distribution
                    self.model_af.loc[:, n] = [item[0] for item in np.random.beta(100 * N_A / N_T, 100 * N_R / N_T)]            
        else:  # final round, initialise with top20 confidently assigned barcodes
            for n in range(1, self.singlets + 1):
                topN = (self.P_s_c[n].argsort() >= (len(self.barcodes) - 20)) * 1
                N_A = self.alt_bc_mtx * topN + self.pseudo
                N_R = self.ref_bc_mtx * topN + self.pseudo
                N_T = N_A + N_R
                k_ref = N_R / N_T
                k_alt = N_A / N_T
                self.model_af.loc[:, n] = (N_A + k_alt) / (N_T + k_alt + k_ref)
                
        ### initialise doublets:
        
        index = self.singlets + 1  # start from doublet states
        for i in range(1, self.singlets):
            for j in range(i + 1, self.singlets + 1):
                self.P_s.append(dbl / (self.num - self.singlets - 1))     # even distribution of doublet probability on all doublet states
                self.model_af.loc[:, index] = (self.model_af.loc[:, i] + self.model_af.loc[:, j]) / 2   # mean of relevant two singlet states
                index += 1


    def calculate_cell_likelihood(self):
        """
        Calculate cell|sample likelihood P(c|s) and derive sample|cell probability P(s|c)
        P(c|s_v) = P(N(A),N(R)|s) = P(g_A|s)^N(A) * (1-P(g_A|s))^N(R)
        log(P(c|s)) = sum_v{(N(A)_c,v*log(P(g_A|s)) + N(R)_c,v*log(1-P(g_A|s)))}
        P(s_n|c) = P(c|s_n) / [P(c|s_1) + P(c|s_2) + ... + P(c|s_n)]
        log(P(s1|c) = log{1/[1+P(c|s2)/P(c|s1)]} = -log[1+P(c|s2)/P(c|s1)] = -log[1+2^(logP(c|s2)-logP(c|s1))]

        """

        # calculate likelihood P(c|s) based on allele probability
        for n in range(self.num):
            matcalc = self.alt_bc_mtx.T.multiply(self.model_af.loc[:, n].apply(np.log2)).T \
                    + self.ref_bc_mtx.T.multiply((1 - self.model_af.loc[:, n]).apply(np.log2)).T
            self.lP_c_s.loc[:, n] = matcalc.sum(axis=0).tolist()[0]  # log likelihood to avoid python computation limit of 1e-323/1e+308
    
        # transform to cell sample probability using Baysian rule
        for i in range(self.num):
            denom = 0
            for j in range(self.num):
                denom += 2 ** (self.lP_c_s.loc[:, j] + np.log2(self.P_s[j]) - self.lP_c_s.loc[:, i] - np.log2(self.P_s[i]))
            self.P_s_c.loc[:, i] = 1 / denom


    def calculate_model_af(self):
        """
        Update the model allele fraction by distributing the alt and total counts of each barcode on a certain snv to the model based on P(s|c)

        """

        N_ref = self.ref_bc_mtx.sum(axis=1) + self.pseudo
        N_alt = self.alt_bc_mtx.sum(axis=1) + self.pseudo
        k_ref = N_ref / (N_ref + N_alt)
        k_alt = N_alt / (N_ref + N_alt)
        self.model_af = pd.DataFrame((self.alt_bc_mtx.dot(self.P_s_c) + k_alt) / ((self.alt_bc_mtx + self.ref_bc_mtx).dot(self.P_s_c) + k_ref + k_alt),
                                        index = self.all_POS, columns = range(self.num))
        self.model_af.loc[:, 0] = self.model_af.loc[:, 1:(self.num-1)].mean(axis=1)   # reset the multiplet AF
        # reset doublet AF
        index = self.singlets + 1  # start from doublet states
        for i in range(1, self.singlets):
            for j in range(i + 1, self.singlets + 1):
                self.model_af.loc[:, index] = (self.model_af.loc[:, i] + self.model_af.loc[:, j]) / 2   # mean of relevant two singlet states
                index += 1

        # reset sample prior probabilities based on sum(P(s|c))
        #self.P_s = self.P_s_c.sum(axis=0).tolist()
        #self.P_s = [item/sum(self.P_s) for item in self.P_s]


    def next_seed(self, x):
        if x <= self.singlets:
            for i in range(2,len(self.barcodes)):
                j = (self.ref_bc_mtx + self.alt_bc_mtx).sum(axis=0).argsort()[0,len(self.barcodes)-i]
                if (max(self.P_s_c.iloc[j, range(1,x+1)]) < 0.5) & (not(j in self.seeds)):
                    self.seeds.append(j)
                    break


    def assign_cells(self):
        """
	    Final assignment of cells according to P(s|c) >= 0.9

	    """

        for n in range(self.num):
            self.assigned[n] = sorted(self.P_s_c.loc[self.P_s_c[n] >= 0.9].index.values.tolist())



def run_model(base_calls_mtx, num_models):

    model = models(base_calls_mtx, num_models)
    
    for m in range(1, num_models + 1):  # num_model rounds, with each round one more additional barcode picked as seed for intialisation

        # initialise model
        model.initialise_model(m)

        # commencing E-M
        iterations = 0
        sum_log_likelihood = [1,2]  # dummy likelihood as a start
        while sum_log_likelihood[-2] != sum_log_likelihood[-1]:
            iterations += 1
            progress = 'Iteration ' + str(iterations) + '   ' + str(datetime.datetime.now()) + '\n'
            with open('wip.log', 'a') as myfile: myfile.write(progress)
            model.calculate_cell_likelihood()  # E-step, calculate the expected cell origin likelihood with a function of model.model_af (theta)
            model.calculate_model_af()  # M-step, to optimise unknown model parameter model.model_af (theta)
            # approximation due to python calculation limit
            sum_log_likelihood.append(model.lP_c_s.max(axis=1).sum())  # L = Prod_c[Sum_s(P(c|s))], thus LL = Sum_c{log[Sum_s(P(c|s))]}
            # sum_log_likelihood.append(((2**model.lP_c_s).sum(axis=1)+1e-323).apply(np.log2).sum())

        model.next_seed(m)

    model.assign_cells()

    # generate outputs
    for n in range(model.num):
        with open('barcodes_{}.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')
    model.P_s_c.to_csv('P_s_c.csv')
    model.model_af.to_csv('model_af.csv')
    for index, item in enumerate(model.P_s_c.index[model.seeds]):
        with open ('seed_{}.txt'.format(str(index+1)), 'w') as myfile:
            myfile.write(item)
    print(model.P_s)
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

    num_models = 4          # number of models in each run

    # input and output files
    ref_csv = 'ref_filtered.csv'  # reference matrix
    alt_csv = 'alt_filtered.csv'  # alternative matrix

    progress = 'Starting data collection: ' + str(datetime.datetime.now()) + '\n'
    with open('wip.log', 'a') as myfile: myfile.write(progress)
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, num_models)

if __name__ == '__main__':
    main()
