"""
Reference free AF-based demultiplexing on pooled scRNA-seq
Genotype verification comparing sc-Split results with known genotypes
Jon Xu (jun.xu@uq.edu.au)
Lachlan Coin
Aug 2018
"""

import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import numpy as np
import pandas as pd
from scipy.stats import binom
from scipy.sparse import csr_matrix
import datetime
import csv


def main():

    # input and output files
    ori_vcf = 'imputed.vcf'        # original vcf file used in demuxlet run
    ref_csv = 'ref_filtered.csv'  # reference matrix
    alt_csv = 'alt_filtered.csv'  # alternative matrix
    psc_csv = 'P_s_c.csv'         # P(s|c) generated by sc_split
    dem_grp = 'd'                 # prefix of demuxlet assignments

    ref = pd.read_csv(ref_csv, header=0, index_col=0)
    alt = pd.read_csv(alt_csv, header=0, index_col=0)
    ref_s = csr_matrix(ref.values)
    alt_s = csr_matrix(alt.values)
    all_POS = ref.index

    # build matrix from scSplit assignment
    P_s_c = pd.read_csv(psc_csv, header=0, index_col=0)
    A_s_c = ((P_s_c >= 0.9) * 1).astype('float64')
    num = len(P_s_c.columns)  # number of samples + 1 doublet state

    err = 0.01  # error rate assumption

    # binomial simulation for genotype likelihoods P(D|AA,RA,RR) with the alt count vs total count condition and (err, 0.5, 1-err) as allele probability
    # cells added to assigned states
    # rr/ra/aa: GL, log likelihood of genotypes based on REF/ALT alleles
    rr = pd.DataFrame(binom.logpmf(pd.DataFrame(alt_s.dot(A_s_c)), pd.DataFrame((alt_s + ref_s).dot(A_s_c)), err), index=all_POS, columns=range(num)).drop(0,1)
    ra = pd.DataFrame(binom.logpmf(pd.DataFrame(alt_s.dot(A_s_c)), pd.DataFrame((alt_s + ref_s).dot(A_s_c)), 0.5), index=all_POS, columns=range(num)).drop(0,1)
    aa = pd.DataFrame(binom.logpmf(pd.DataFrame(alt_s.dot(A_s_c)), pd.DataFrame((alt_s + ref_s).dot(A_s_c)), 1-err), index=all_POS, columns=range(num)).drop(0,1)

    # prr/pra/paa: GP, genotype posterior probabilities P(AA,RA,RR|D)
    # transform to cell sample probability using Baysian rule
    prr = 1 / (1 + (ra - rr).apply(np.exp) + (aa - rr).apply(np.exp))
    pra = 1 / (1 + (rr - ra).apply(np.exp) + (aa - ra).apply(np.exp))
    paa = 1 / (1 + (rr - aa).apply(np.exp) + (ra - aa).apply(np.exp))

    # locate most likely genotypes at each SNV for each sample and build genotype matrix for comparison
    prr = (prr > 0.99) * 1      # start with 1 to avoid 0 confusion
    pra = (pra > 0.99) * 2
    paa = (paa > 0.99) * 3
    scsplit = prr + pra + paa   # for each SNV and sample, the value will only be 1, 2, 3 or NA

    # build genotype matrix from vcf                                    ### are SNVs in the ori_vcf same as the ref/alt matrices
    prr[:] = pra[:] = paa[:] = 0    
    for record in vcf.Reader(open(ori_vcf, 'r')):
        pos = record.CHROM + ':' + str(record.POS)
        if pos in list(prr.index):
            for n in range (1, num):
                if record.samples[n-1]['GP'][0] > 0.99:
                    prr.loc[pos, n] = 1
                elif record.samples[n-1]['GP'][1] > 0.99:
                    prr.loc[pos, n] = 2
                elif record.samples[n-1]['GP'][2] > 0.99:
                    prr.loc[pos, n] = 3
    vcf_all = prr + pra + paa

    # output
    scsplit.to_csv('verify_scsplit.csv')
    vcf_all.to_csv('verify_vcf.csv')

if __name__ == '__main__':
    main()
