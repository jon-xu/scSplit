"""
Build SNV/Barcode Matrices (ref/alt) for genotype-free demultiplexing on pooled scRNA-seq
Jon Xu
Lachlan Coin
Aug 2018
jun.xu@uq.edu.au
"""

import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import math, argparse, datetime
import numpy as np
import pysam as ps  # http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
import pandas as pd

class SNV_data:
    """
    Stores data on each SNV
    """

    def __init__(self, chrom, pos, ref, alt):

        self.CHROM = chrom   # chromosome number
        self.POS = pos       # position on chromosome
        self.REF = ref       # reference base
        self.ALT = alt       # alternate base


def build_base_calls_matrix(file_s, all_SNVs, barcodes):
    """
    Build pandas DataFrame
    Parameters:
        file_s(str): Path to sam file (0-based positions)
        all_SNVs: list of SNV_data objects
        barcodes(list): cell barcodes
    """

    all_POS = []   # snv positions (1-based positions from vcf file)
    for entry in all_SNVs:
        pos = str(entry.CHROM) + ':' + str(entry.POS)
        if pos not in all_POS:
            all_POS.append(pos)

    in_sam = ps.AlignmentFile(file_s, 'rb')
    ref_base_calls_mtx = pd.DataFrame(0, index=all_POS, columns=barcodes)
    alt_base_calls_mtx = pd.DataFrame(0, index=all_POS, columns=barcodes)
    print('Num Pos:', len(all_POS), ', Num barcodes:', len(barcodes))

    for snv in all_SNVs:
        position = str(snv.CHROM) + ':' + str(snv.POS)
        # use pysam.AlignedSegment.fetch instead of pysam.AlignedSegment.pileup which doesn't contain barcode information
        for read in in_sam.fetch(snv.CHROM, snv.POS-1, snv.POS+1):
            if read.flag < 256:   # only valid reads
                if (snv.POS - 1) in read.get_reference_positions():
                    # if the read aligned positions cover the SNV position
                    try:
                        barcode = read.get_tag('CB')
                    except:
                        barcode = ''
                    if barcode in barcodes:
                        # read the base from the snv.POS which the read has mapped to
                        base = read.query_sequence[[item for item in read.get_aligned_pairs(True) if item[1] == (snv.POS - 1)][0][0]]
                        if base == snv.REF:
                            ref_base_calls_mtx.loc[position, barcode] += 1
                        if base == snv.ALT:
                            alt_base_calls_mtx.loc[position, barcode] += 1

    ref_base_calls_mtx.index.name = alt_base_calls_mtx.index.name = 'SNV'

    return (ref_base_calls_mtx, alt_base_calls_mtx)


def main():

    # Process command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', required=True,  help='VCF from mixed BAM')
    parser.add_argument('-i', '--bam', required=True, help='mixed sample BAM')
    parser.add_argument('-b', '--barcodes', required=True,  help='barcodes')
    parser.add_argument('-r', '--ref', required=True,  help='Ref count CSV')
    parser.add_argument('-a', '--alt', required=True,  help='Alt count CSV')
 
    args = parser.parse_args()
    dist_alleles = []
    
    epsilon = 0.01

    all_SNVs = []  # list of SNV_data objects
    for record in vcf.Reader(open(args.vcf, 'r')):
        # only keep SNVs with heterozygous genotypes, and ignore SNV with multiple bases (e.g. GGGT/GGAT)
        if (record.samples[0]['GL'][1] > np.log10(1-epsilon)) & (len(record.REF) == 1) & (len(record.ALT) == 1):
            all_SNVs.append(SNV_data(record.CHROM, record.POS, record.REF, record.ALT[0]))
    
    barcodes = []   # list of cell barcodes
    for line in open(args.barcodes, 'r'):
        barcodes.append(line.strip())

    base_calls_mtx = build_base_calls_matrix(args.bam, all_SNVs, barcodes)
    base_calls_mtx[0].to_csv('{}'.format(args.ref))
    base_calls_mtx[1].to_csv('{}'.format(args.alt))

if __name__ == '__main__':
    main()
