"""
Reference free genotype demultiplexing on pooled scRNA-seq
Jon Xu, Caitlin Falconer
Aug 2018
"""

import pdb
import sys
import numpy as np
import pandas as pd
import datetime
import csv
from scipy.stats import binom

BASE_RR_PROB = BASE_RA_PROB = BASE_AA_PROB = 0.3333 

class SNV_data:
    """
    Stores data on each SNV
    """

    def __init__(self, chrom_pos):
        """
        Parameters:
            chrom (int): chromosome number
            pos (int): position on chromosome
        """
        self.CHROM = chrom_pos.split(':')[0] 
        self.POS = chrom_pos.split(':')[1]

    def get_all_SNV_pos(all_SNVs):
        """
        Return ordered list of positions of SNVs as chr:pos
        Parameters:
            all_SNVs (list(SNV_data obj)): list of SNV_data objects
        Returns:
            sorted list of SNV unique positions as chr:pos
        """
        all_POS = []
        for entry in all_SNVs:
            pos = str(entry.CHROM) + ':' + str(entry.POS)
            if pos not in all_POS:
                all_POS.append(pos)
        return all_POS

class models:
    def __init__(self, all_SNVs, base_calls_mtx, barcodes, num=2, model_genotypes=[], assigned=None, P_c_s=[], P_s_c=[]):
        """
        A complete model class containing SNVs, matrices counts, barcodes, model genotypes with assigned cells 

        Parameters:
             all_SNVs: list[SNV_data objects]
             base_calls_mtx: SNV-barcode matrix containing lists of base calls
             num(int): number of individual model genotypes
             barcodes: list of cell barcodes
             model_genotypes(list(DataFrame): list of num model genotypes represented in snv-probdistrib DataFrame as <RR,RA,AA>
             P_s_c: barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
             P_c_s: barcode/sample matrix containing the probability of seeing barcode c under sample s, whose sum should increase by iterations, as a model indicator 
             assigned(list): final lists of cell/barcode assigned to each genotype cluster/model	     
        """
        self.all_SNVs = all_SNVs
        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.barcodes = barcodes
        self.num = num
        self.model_genotypes = model_genotypes
        self.assigned = assigned
        self.P_s_c = P_s_c
        self.P_c_s = P_c_s

        if self.P_s_c == []:
            self.P_s_c = pd.DataFrame(np.zeros((len(self.barcodes), num)),
                            index = self.barcodes, columns = list(range(num)))

        if self.P_c_s == []:
            self.P_c_s = pd.DataFrame(np.zeros((len(self.barcodes), num)),
                            index = self.barcodes, columns = list(range(num)))

        if self.assigned == None:
            self.assigned = []
            for _ in range(self.num):
                self.assigned.append([])


    def calculate_model_genotypes(self):
        """
	    Initialise model from dirichlet distribution
	    """
        if self.model_genotypes == []:
            for n in range(self.num):
                self.model_genotypes.append([])
                self.model_genotypes[n] = pd.DataFrame(np.zeros((len(self.all_SNVs), 3)),
                    index=SNV_data.get_all_SNV_pos(self.all_SNVs), columns=['RR', 'RA', 'AA'])

            for n in range(self.num):
                dir_sim = np.random.dirichlet((1,1),len(self.all_SNVs))  # generate random dirichlet distribution to simulate genotypes probability
                gt = [item[0] for item in dir_sim]   # take the first value of each pair in the distribution results
                i = 0
                for snv in self.all_SNVs:
                    pos = "{}:{}".format(snv.CHROM, snv.POS)
                    self.model_genotypes[n].loc[pos] = ((1-gt[i])**2, 2*gt[i]*(1-gt[i]), gt[i]**2)   # based on Hardy-Weinberg equilibrium
                    i += 1

        else:
            """
            Update the probability distribution <RR, RA, AA> for each SNV position based on counts of bases

            Input:
                bc_mtx: SNV-barcode matrix containing lists of base calls
                P_s_c: the cell likelihood from last iteration, weight the alt count distribution in calculating the genotype

            Output:
                model_genotypes: P(g|S_v) according to the alt counts and the P(s|c)
            """

            A_sv = T_sv = 0     # P(s|c) weighted Alt/Total matrices counts for the SNV in sample
            for snv in self.all_SNVs:  # for each snv position
                pos = "{}:{}".format(snv.CHROM, snv.POS)
                for n in range(self.num):  # for the number of genotypes being modelled
                    for barcode in self.barcodes:
                        A_sv += self.P_s_c.loc[barcode, n] * self.alt_bc_mtx.loc[pos, barcode]
                        T_sv += self.P_s_c.loc[barcode, n] * (self.alt_bc_mtx.loc[pos, barcode] + self.ref_bc_mtx.loc[pos, barcode])
                    
                    if A_sv > 0 or T_sv > 0:      # if cells assigned in model carry reads at this snv
                        P_A_S = (A_sv + 1) / (T_sv + 2)       # probability of alt counts, with pseudo counts
                        self.model_genotypes[n].loc[pos] = ((1-P_A_S)**2, 2*P_A_S*(1-P_A_S), P_A_S**2)   # derive genotype probabilities using observed ALT probability, HWE


    def calculate_cell_likelihood(self):
        """
        Allocate alternative counts to different genotype clusters within a model
        Calculate likelihood cell/barcode came from individual i (i.e. model genotype 1...n):
        For each cell/barcode c, product over all variants v, sum over all genotypes <RR,RA,AA>
        Product over all reads in cell c for variant v
        Pr(base_called | g)*P(g|s)
        where P(g|s) = set(RR,RA,AA)
        
        Input:
            genotype model: self.model_genotypes
            snv pos: self.model_genotypes[0].index
            bc_mtx: SNV-barcode matrix containing lists of base calls
        
        Output:
            P_s_c: barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
            P_c_s: barcode/sample matrix containing the probability of seeing barcode c under sample s, whose sum should increase by iterations, as a model indicator 

        Issue: 
            Some barcodes in matrix have no base calls as they do not cover snp region,
            due to base_calls_mtx being instantiated from list of all barcodes in bam file (not just those covering snv)
        """

        err = 0.01

        for barcode in self.barcodes:  # for each of all cells
            cell_llhood = []  # store likelihood values of cell/barcode c for model 1 to n

            indices = []  # index of non-zero values in column
            for idx in self.ref_bc_mtx.loc[:, barcode].nonzero():
                indices.extend(idx)
            for idx in self.alt_bc_mtx.loc[:, barcode].nonzero():
                indices.extend(idx)
            for n in range(self.num):  # for each model
                all_var_llhood = 1  # initialise the product over all variants for model n
                for i in indices:  # for each SNV with non-zero values in this cell
                    snv = self.all_SNVs[i]
                    pos = "{}:{}".format(snv.CHROM, snv.POS)

                    ref_count = self.ref_bc_mtx.loc[pos, barcode]   # count of ref reads at a SNV position from a cell
                    alt_count = self.alt_bc_mtx.loc[pos, barcode]
                    tot_count = ref_count + alt_count

                    # calculate probability of alt_count in the barcode at the SNV conditioned on general genotype
                    # probability of A in RR/RA/AA is 0.01, 0.5, 0.99
                    P_A_g = binom.pmf(alt_count, tot_count, (err, 0.5, 1-err))

                    # calculate P(c|S_n) for the SNV, sum of the three P_c_given_g,S for model n
                    P_c_Sv = sum(P_A_g * self.model_genotypes[n].loc[pos])

                    # Product over all variant positions
                    all_var_llhood *= P_c_Sv

                self.P_c_s.loc[barcode,n] = all_var_llhood  # array of [P(c|S_1), P(c|S_2)] for each model

            # transform the P(c|S_n) of all cells into posterior probability, assuming P(S_n) are all equal (prior probability)
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


def run_model(all_SNVs, base_calls_mtx, barcodes, num_models):

    model = models(all_SNVs, base_calls_mtx, barcodes, num_models, model_genotypes=[], assigned=None, P_c_s=[], P_s_c=[])

    print("Commencing E-M")
    
    iterations = 0
    sum_log_likelihoods = []

    while iterations < 20:

        iterations += 1

        print("calculating model ", datetime.datetime.now().time())
        model.calculate_model_genotypes()

        print("calculating cell likelihood ", datetime.datetime.now().time())
        model.calculate_cell_likelihood()

        sum_log_likelihood = model.P_c_s.apply(np.log).sum().sum()
        sum_log_likelihoods.append(sum_log_likelihood)
        print("log likelihood of iteration {}".format(iterations), sum_log_likelihood)

    model.assign_cells()

    for n in range(num_models):
        with open('barcodes_{}_de.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')

        model.model_genotypes[n].to_csv('model_genotypes{}_de.csv'.format(n))

    model.P_c_s.to_csv('P_c_s.csv')
    model.P_s_c.to_csv('P_s_c.csv')
    
    print("Finished model at {}".format(datetime.datetime.now().time()))
    print(sum_log_likelihoods)

def get_all_barcodes(bc_file):
    """
    Load barcodes from text file
    Parameters:
        bc_file: text file of line separated cell barcodes
    Returns:
        list of cell barcodes
    """
    file = open(bc_file, 'r')
    barcodes = []
    for line in file:
        barcodes.append(line.strip())
    return barcodes


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
    #bc_file = "bc_sorted.txt"   # validated cell barcodes
    #ref_csv = 'ac21_ref.csv'  # reference matrix
    #alt_csv = 'ac21_alt.csv'  # alternative matrix
    #outfile_A = 'model_A{}_21.txt'
    #outfile_B = 'model_B{}_21.txt' 

    bc_file = 'test.txt'
    ref_csv = 'test_ref.csv'
    alt_csv = 'test_alt.csv'
    outfile_A = 'model_A{}_de.txt'
    outfile_B = 'model_B{}_de.txt'
    outfile_g = 'model_genotypes{}{}_de.csv'

    all_SNVs = []  # list of SNV_data objects
    matrix = pd.read_csv(ref_csv)
    for record in matrix.loc[:,'SNV']:
        all_SNVs.append(SNV_data(record))
    
    print ("Starting data collection", datetime.datetime.now().time())
    
    all_POS = SNV_data.get_all_SNV_pos(all_SNVs)
    barcodes = get_all_barcodes(bc_file)
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(all_SNVs, base_calls_mtx, barcodes, num_models)

if __name__ == "__main__":
    main()
