#!/usr/bin/env python3

"""
Reference free AF-based demultiplexing on pooled scRNA-seq
Genotype validation comparing sc-Split results with known genotypes
Jon Xu (jun.xu@uq.edu.au)
Lachlan Coin
Aug 2018
"""

import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import numpy as np
import pandas as pd
from scipy.stats import binom
from scipy.sparse import csr_matrix
import argparse, datetime, csv


def main():

    # Process command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', required=True,  help='Ref count CSV')
    parser.add_argument('-a', '--alt', required=True,  help='Alt count CSV')
    parser.add_argument('-i', '--input', required=True,  help='distinguishing alleles')
    parser.add_argument('-v', '--vcf', required=True,  help='imputed VCF')
    parser.add_argument('-p', '--psc', required=True, help='generated P(S|C)')
    parser.add_argument('-m', '--matrix', required=True,  help='imputed matrix')
    args = parser.parse_args()
    dist_alleles = []

    ref = pd.read_csv(args.ref, header=0, index_col=0)
    alt = pd.read_csv(args.alt, header=0, index_col=0)
    ref_s, alt_s = csr_matrix(ref.values), csr_matrix(alt.values)
    all_POS = ref.index    

    ### check genotypes of distinguishing alleles according to sample VCF
    for line in open(args.input, 'r'):
        dist_alleles.append(line.strip())
    matrix_imputed = pd.DataFrame(-1, index=dist_alleles, columns=range(8))
    for item in dist_alleles:
        for record in vcf.Reader(open(args.vcf, 'r')):
            if str(record.CHROM) == item.split(':')[0]:
                if str(record.POS) == item.split(':')[1]:
                    for s, sample in enumerate(record.samples):
                        if sample['GP'][0] > 0.99:
                            matrix_imputed.loc[item][s] = 0
                        elif sample['GP'][0] < 0.01:
                            matrix_imputed.loc[item][s] = 1
                    break

    matrix_imputed.to_csv(args.matrix)

    ### check genotype similarity

    # build matrix from scSplit assignment
    P_s_c = pd.read_csv(args.psc, header=0, index_col=0)
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

    # build genotype matrix from vcf
    prr[:] = pra[:] = paa[:] = 0    
    for record in vcf.Reader(open(args.vcf, 'r')):
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
    scsplit.to_csv('scsplit_validate_scsplit.csv')
    vcf_all.to_csv('scsplit_validate_vcf.csv')

if __name__ == '__main__':
    main()
