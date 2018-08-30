"""
Build SNV/Barcode Matrices (ref/alt) for genotype demultiplexing on pooled scRNA-seq
"""

import pdb
import sys
import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import numpy as np
import pysam as ps  # http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
import pandas as pd
import subprocess
import datetime

class SNV_data:
    """
    Stores data on each SNV
    """

    def __init__(self, chrom, pos, ref, alt):
        """
        Parameters:
            chrom (int): chromosome number
            pos (int): position on chromosome
            ref (str): reference base
            alt (str): alternate base
            loglik (tuple(float)): log likelihood of genotypes normalised to most likely allele (RR,RA,AA)
        """
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alt

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


def build_base_calls_matrix(sam_filename, all_SNVs, all_POS, barcodes):
    """
    Build pandas DataFrame
    Parameters:
        sam_filename(str): Path to sam file (0-based positions)
        all_SNVs: list of SNV_data objects
        all_POS(list('chr:pos)): snv positions (1-based positions from vcf file)
        barcodes(list): cell barcodes
    """

    in_sam = ps.AlignmentFile(sam_filename, 'rb')
    ref_base_calls_mtx = pd.DataFrame(np.zeros((len(all_POS), len(barcodes))),
                                      index=all_POS, columns=barcodes)
    alt_base_calls_mtx = pd.DataFrame(np.zeros((len(all_POS), len(barcodes))),
                                      index=all_POS, columns=barcodes)
    print('Matrix size: Num Pos:', len(all_POS), 'Num barcodes:', len(barcodes))

    all_POS = []
    for entry in all_SNVs:
        pos = str(entry.CHROM) + ':' + str(entry.POS)
        if pos not in all_POS:
            all_POS.append(pos)
    for snv in all_SNVs:
        position = str(snv.CHROM) + ':' + str(snv.POS)
        # use pysam.AlignedSegment.fetch to replace pysam.AlignedSegment.pileup which doesn't contain barco    de information
        for read in in_sam.fetch(snv.CHROM, snv.POS-1, snv.POS+1):
            if read.flag < 256:   # only valid reads
                if (snv.POS - 1) in read.get_reference_positions():
                     # if the read aligned positions cover the SNV position
                    try:
                        barcode = read.get_tag('CB')
                    except:
                        barcode = ''
                    if barcode is not '':
                        # read the base from the snv.POS which the read has mapped to
                        base = read.query_sequence[[item for item in read.get_aligned_pairs(True) if item[1] == (snv.POS - 1)][0][0]]
                        if base == snv.REF:
                            ref_base_calls_mtx.loc[position, barcode] += 1
                        if base == snv.ALT:
                            alt_base_calls_mtx.loc[position, barcode] += 1

    ref_base_calls_mtx.index.name = alt_base_calls_mtx.index.name = 'SNV'

    return (ref_base_calls_mtx, alt_base_calls_mtx)


def main():

    """ The following is various files used while testing the program
    file_v = .vcf file
    file_s = .bam file
    file_bc = .txt file of all known and checked barcodes
    """

    # Mixed donor files
    file_v = "ac21.vcf"
    file_s = "ac21_q10_rmdup_sorted.bam"
    file_bc = "bc_sorted.txt"
    out_csv_ref = 'ac21_ref.csv'
    out_csv_alt = 'ac21_alt.csv'

    in_vcf = vcf.Reader(open(file_v, 'r'))

    vcf_records = []
    for record in in_vcf:
        vcf_records.append(record)

    all_SNVs = []  # list of SNV_data objects
    for record in vcf_records:
        if record.samples[0]['GL'][0] != 0 and record.samples[0]['GL'][2] != 0:
            all_SNVs.append(
                SNV_data(record.CHROM, record.POS, record.REF, record.ALT[0]))

    print("Starting data collection", datetime.datetime.now().time())
    all_POS = SNV_data.get_all_SNV_pos(all_SNVs)
    barcodes = get_all_barcodes(file_bc)

    base_calls_mtx = build_base_calls_matrix(file_s, all_SNVs,
                                             all_POS,barcodes)
    base_calls_mtx[0].to_csv('{}{}{}'.format(path, out_dir, out_csv_ref))
    base_calls_mtx[1].to_csv('{}{}{}'.format(path, out_dir, out_csv_alt))
    print("Base call matrix finished", datetime.datetime.now().time())


if __name__ == "__main__":
    main()
