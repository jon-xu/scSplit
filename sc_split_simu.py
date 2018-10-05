"""
Simulate SNV/Barcode Matrices (ref/alt) for genotype-free demultiplexing on pooled scRNA-seq Jon Xu Lachlan Coin Aug 2018 jun.xu@uq.edu.au """

import sys
import vcf  # https://pyvcf.readthedocs.io/en/latest/INTRO.html
import random
import numpy as np
import pysam as ps  # http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
import pandas as pd
import datetime

class SNV_data:
    """
    Stores data on each SNV
    """

    def __init__(self, chrom, pos, ref, alt, samples):

        self.CHROM = chrom   # chromosome number
        self.POS = pos       # position on chromosome
        self.REF = ref       # reference base
        self.ALT = alt       # alternate base
        self.SAMPLES = samples   # mixed samples


def simulate_base_calls_matrix(file_i, file_o, all_SNVs, barcodes):
    """
    Build pandas DataFrame
    Parameters:
        file_s(str): Path to sam file (0-based positions)
        all_SNVs: list of SNV_data objects
        barcodes(list): cell barcodes
    """

    num = 2  # user input

    # randomly put all barcodes into the groups
    groups = [random.randint(0, num-1) for item in range(len(barcodes))]

    # output randomly assigned cell groups
    with open('bc_groups.txt', 'w+') as myfile:
        for n in range(len(barcodes)):
            myfile.write(barcodes[n] + '\t' + str(groups[n]) + '\n')

    all_POS = []   # snv positions (1-based positions from vcf file)
    for entry in all_SNVs:
        pos = str(entry.CHROM) + ':' + str(entry.POS)
        if pos not in all_POS:
            all_POS.append(pos)

    in_sam = ps.AlignmentFile(file_i, 'rb')
    out_sam = ps.AlignmentFile(file_o, 'wb', template=in_sam)

    ref_base_calls_mtx = pd.DataFrame(np.zeros((len(all_POS), len(barcodes))), index=all_POS, columns=barcodes)
    alt_base_calls_mtx = pd.DataFrame(np.zeros((len(all_POS), len(barcodes))), index=all_POS, columns=barcodes)
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

                        # get sample id from randomly generated grouping indexes for the list of barcodes
                        sample = groups[barcodes.index(barcode)]

                        # calculate probability of A allele based on the genotype probability (GP) or likelihoods (10 ^ GL)
                        try:
                            rr = 10 ** snv.SAMPLES[sample]['GL'][0]
                            ra = 10 ** snv.SAMPLES[sample]['GL'][1]
                            aa = 10 ** snv.SAMPLES[sample]['GL'][2]
                            # transform P(D|G) to P(G|D)
                            P_A = (ra/2 + aa) / (rr + ra + aa)
                            # P(A) = 0.5 P(RA) + P(AA)
                        except:
                            P_A = 0.5 * snv.SAMPLES[sample]['GP'][1] + snv.SAMPLES[sample]['GP'][2]

                        # toss a biased coin using P_A to get A/R allele for the simulated read
                        if random.random() < P_A:
                            alt_base_calls_mtx.loc[position, barcode] += 1
                            new = str(snv.ALT[0])  # ALT returned as list by pysam
                        else:
                            ref_base_calls_mtx.loc[position, barcode] += 1
                            new = snv.REF

                        # update the new base to bam file
                        for item in read.get_aligned_pairs(True):
                            if item[1] == (snv.POS - 1):
                                read.query_sequence = read.query_sequence[:item[0]] + new + read.query_sequence[(item[0] + 1):]
                        out_sam.write(read)

                        # simulate doublets
                        if barcodes.index(barcode) < len(barcodes) * 0.02:  # assume 2% doublets
                            sample = (groups[barcodes.index(barcode)] + 1 ) % num
                            P_A = 0.5 * snv.SAMPLES[sample]['GP'][1] + snv.SAMPLES[sample]['GP'][2]
                            if random.random() < P_A:
                                alt_base_calls_mtx.loc[position, barcode] += 1
                                new = str(snv.ALT[0])
                            else:
                                ref_base_calls_mtx.loc[position, barcode] += 1
                                new = snv.REF
                            for item in read.get_aligned_pairs(True):
                                if item[1] == (snv.POS - 1):
                                    read.query_sequence = read.query_sequence[:item[0]] + new + read.query_sequence[(item[0] + 1):]
                            out_sam.write(read)

                    except:
                        pass
                
                else:
                    out_sam.write(read)  # write the read to out_sam without changes if SNV position is not a reference position
            
            else:
                out_sam.write(read)  # write the read to out_sam without changes if its flag >= 256

        # reads not overlapping with SNV position is filtered out

    ref_base_calls_mtx.index.name = alt_base_calls_mtx.index.name = 'SNV'

    return (ref_base_calls_mtx, alt_base_calls_mtx)


def main():

    # Input and output files
    file_v = 'input.vcf'
    file_b = 'input.txt'
    file_i = 'input.bam'
    file_o = 'sim.bam'
    out_csv_ref = 'ref_sim.csv'
    out_csv_alt = 'alt_sim.csv'
    
    all_SNVs = []  # list of SNV_data objects
    for record in vcf.Reader(open(file_v, 'r')):
        # haven't decided which sample here, just ignore SNV with multiple bases (e.g. GGGT/GGAT)
        if (len(record.REF) == 1) & (len(record.ALT) == 1):
            all_SNVs.append(SNV_data(record.CHROM, record.POS, record.REF, record.ALT, record.samples))
    
    barcodes = []   # list of cell barcodes
    for line in open(file_b, 'r'):
        barcodes.append(line.strip())

    base_calls_mtx = simulate_base_calls_matrix(file_i, file_o, all_SNVs, barcodes)
    base_calls_mtx[0].to_csv('{}'.format(out_csv_ref))
    base_calls_mtx[1].to_csv('{}'.format(out_csv_alt))

if __name__ == "__main__":
    main()
