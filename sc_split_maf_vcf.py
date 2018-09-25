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
        self.model_genotypes = []
        for _ in range(self.num):
            self.assigned.append([])
            self.model_genotypes.append([])
        self.model_MAF = pd.DataFrame(np.zeros((len(self.all_POS), self.num)), index=self.all_POS, columns=range(self.num))
        # set background alt count proportion as fixed minor allele frequency for each SNVs in the model
        self.model_MAF.loc[:, 0] = self.alt_bc_mtx.sum(axis=1) / (self.ref_bc_mtx.sum(axis=1) + self.alt_bc_mtx.sum(axis=1))
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
            self.lP_c_s.loc[:, n] = matcalc.sum(axis=0).tolist()[0]  # log likelihood to avoid python computation limit of 2^+/-308
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


    def calculate_model_genotypes(self):
        """
        Update the model genotype by simulating the count distribution using P(s|c) and P(g|D) of assigned barcode on a certain snv to the final model

        """

        err = 0.01  # error rate assumption
        # binomial probability for P(D|AA,RA,RR) with the alt count vs total count condition and (err, 0.5, 1-err) genotype probability
        self.p_d_aa = pd.DataFrame(binom.pmf(pd.DataFrame(self.alt_bc_mtx.todense()), pd.DataFrame((self.alt_bc_mtx + self.ref_bc_mtx).todense()), 1-err), index = self.all_POS, columns = self.barcodes)
        self.p_d_ra = pd.DataFrame(binom.pmf(pd.DataFrame(self.alt_bc_mtx.todense()), pd.DataFrame((self.alt_bc_mtx + self.ref_bc_mtx).todense()), 0.5), index = self.all_POS, columns = self.barcodes)
        self.p_d_rr = pd.DataFrame(binom.pmf(pd.DataFrame(self.alt_bc_mtx.todense()), pd.DataFrame((self.alt_bc_mtx + self.ref_bc_mtx).todense()), err), index = self.all_POS, columns = self.barcodes)

        # transform into posterior probability P(AA,RA,RR|D)
        self.p_aa_d = self.p_d_aa / (self.p_d_aa + self.p_d_ra + self.p_d_rr)
        self.p_ra_d = self.p_d_ra / (self.p_d_aa + self.p_d_ra + self.p_d_rr)
        self.p_rr_d = self.p_d_rr / (self.p_d_aa + self.p_d_ra + self.p_d_rr)
        
        self.A_s_c = (self.P_s_c >= 0.9) * 1

        self.model_genotypes[0] = self.p_aa_d.dot(self.A_s_c) / (self.p_aa_d.dot(self.A_s_c) + self.p_ra_d.dot(self.A_s_c) + self.p_rr_d.dot(self.A_s_c) + 1)
        self.model_genotypes[1] = self.p_ra_d.dot(self.A_s_c) / (self.p_aa_d.dot(self.A_s_c) + self.p_ra_d.dot(self.A_s_c) + self.p_rr_d.dot(self.A_s_c) + 1)
        self.model_genotypes[2] = self.p_rr_d.dot(self.A_s_c) / (self.p_aa_d.dot(self.A_s_c) + self.p_ra_d.dot(self.A_s_c) + self.p_rr_d.dot(self.A_s_c) + 1)


def run_model(base_calls_mtx, num_models):

    model = models(base_calls_mtx, num_models)
    
    iterations = 0
    sum_log_likelihood = []

    # commencing E-M
    while iterations < 15:
        iterations += 1
        print("Iteration {}".format(iterations))
        model.calculate_cell_likelihood()  # E-step, calculate the expected cell origin likelihood with a function of model.model_MAF (theta)
        model.calculate_model_MAF()  # M-step, to optimise unknown model parameter model.model_MAF (theta)
        sum_log_likelihood.append((2**model.lP_c_s).sum(axis=1).apply(np.log2).sum())

    model.assign_cells()
    model.calculate_model_genotypes()

    # generate outputs
    for n in range(num_models+1):
        with open('barcodes_maf_{}.csv'.format(n), 'w') as myfile:
            for item in model.assigned[n]:
                myfile.write(str(item) + '\n')
    model.P_s_c.to_csv('P_s_c_maf.csv')
    model.model_MAF.to_csv('model_MAF.csv')
    model.model_genotypes[0].to_csv('model_genotypes_AA.csv')  # Likelihood of AA in all models
    model.model_genotypes[1].to_csv('model_genotypes_RA.csv')  # Likelihood of RA in all models
    model.model_genotypes[2].to_csv('model_genotypes_RR.csv')  # Likelihood of RR in all models
    generate_vcf(model)
    print(sum_log_likelihood)
    print("Finished model at {}".format(datetime.datetime.now().time()))


def read_base_calls_matrix(ref_csv, alt_csv):

    """ Read in existing matrix from the csv files """

    ref = pd.read_csv(ref_csv, header=0, index_col=0)
    alt = pd.read_csv(alt_csv, header=0, index_col=0)
    ref_s = csr_matrix(ref.values)
    alt_s = csr_matrix(alt.values)
    base_calls_mtx = [ref_s, alt_s, ref.index, ref.columns]
    print("Base call matrix finished", datetime.datetime.now().time())
    return base_calls_mtx


def generate_vcf(model):
    num = model.num
    vcf_content = pd.DataFrame(index = model.all_POS, columns = range(-8, num))
    names = vcf_content.columns.tolist()
    names[0] = '#CHROM'
    names[1] = 'POS'
    names[2] = 'ID'
    names[3] = 'REF'
    names[4] = 'ALT'
    names[5] = 'QUAL'
    names[6] = 'FILTER'
    names[7] = 'INFO'
    names[8] = 'FORMAT'
    vcf_content.columns = names
    vcf_content.loc[:,'#CHROM'] = [item.split(':')[0] for item in model.all_POS]
    vcf_content.loc[:,'POS'] = [item.split(':')[1] for item in model.all_POS]
    vcf_content.loc[:,'ID'] = model.all_POS
    vcf_content.loc[:,'REF'] = '.'
    vcf_content.loc[:,'ALT'] = '.'
    vcf_content.loc[:,'QUAL'] = '.'
    vcf_content.loc[:,'FILTER'] = '.'
    vcf_content.loc[:,'INFO'] = '.'
    vcf_content.loc[:,'FORMAT'] = 'GL'

    AA = round(model.model_genotypes[0].astype(float) + 0.0005, 3).astype(str)
    RA = round(model.model_genotypes[1].astype(float) + 0.0005, 3).astype(str)
    RR = round(model.model_genotypes[2].astype(float) + 0.0005, 3).astype(str)

    for n in range(1, num):
        vcf_content.loc[:, n] = RR.loc[:, n] + ',' + RA.loc[:, n] + ',' + AA.loc[:, n]

    with open('sc_split.vcf') as myfile:
        myfile.write('##fileformat=VCFv4.2\n')
        myfile.write('##fileDate=')
        myfile.write(str(datetime.datetime.today()).split(' ')[0])
        myfile.write('\n')
        myfile.write('##source=sc_split\n')
        myfile.write('##reference=hg19.fa\n')
        myfile.write('##contig=<ID=1,length=249250621>\n')
        myfile.write('##contig=<ID=10,length=135534747>\n')
        myfile.write('##contig=<ID=11,length=135006516>\n')
        myfile.write('##contig=<ID=12,length=133851895>\n')
        myfile.write('##contig=<ID=13,length=115169878>\n')
        myfile.write('##contig=<ID=14,length=107349540>\n')
        myfile.write('##contig=<ID=15,length=102531392>\n')
        myfile.write('##contig=<ID=16,length=90354753>\n')
        myfile.write('##contig=<ID=17,length=81195210>\n')
        myfile.write('##contig=<ID=18,length=78077248>\n')
        myfile.write('##contig=<ID=19,length=59128983>\n')
        myfile.write('##contig=<ID=2,length=243199373>\n')
        myfile.write('##contig=<ID=20,length=63025520>\n')
        myfile.write('##contig=<ID=21,length=48129895>\n')
        myfile.write('##contig=<ID=22,length=51304566>\n')
        myfile.write('##contig=<ID=3,length=198022430>\n')
        myfile.write('##contig=<ID=4,length=191154276>\n')
        myfile.write('##contig=<ID=5,length=180915260>\n')
        myfile.write('##contig=<ID=6,length=171115067>\n')
        myfile.write('##contig=<ID=7,length=159138663>\n')
        myfile.write('##contig=<ID=8,length=146364022>\n')
        myfile.write('##contig=<ID=9,length=141213431>\n')
        myfile.write('##contig=<ID=MT,length=16569>\n')
        myfile.write('##contig=<ID=X,length=155270560>\n')
        myfile.write('##contig=<ID=Y,length=59373566>\n')
        myfile.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        myfile.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total Allele Count">\n')
        myfile.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate Allele Count">\n')
        myfile.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Alternate Allele Frequency">\n')
        myfile.write('##FORMAT=<ID=GL,Number=3,Type=Float,Description="Genotype Likelihood for RR/RA/AA">\n')
        vcf_content.to_csv(myfile, index=False, sep='\t', quoting=csv.QUOTE_NONE, escapechar='"')

def main():

    num_models = 2          # number of models in each run

    # file names
    ref_csv = 'ref_filtered.csv'  # reference matrix
    alt_csv = 'alt_filtered.csv'  # alternative matrix

    print ("Starting data collection", datetime.datetime.now().time())
    
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    run_model(base_calls_mtx, num_models)

if __name__ == "__main__":
    main()


