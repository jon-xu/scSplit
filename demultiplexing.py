"""
Reference free genotype demultiplexing on pooled scRNA-seq
New algorithm applied by Jon Xu on the original framework of Caitlin's script
"""
 
import pdb
import sys
import numpy as np
from scipy.stats import binom
import pysam as ps  # http://pysam.readthedocs.io/en/latest/api.html#sam-bam-files
import pandas as pd
import datetime

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


class model_genotype:
    def __init__(self, all_SNVs, base_calls_mtx, barcodes, num=2,
                 model_genotypes=[], assigned=None):
        """
        Parameters:
             all_SNVs: list[SNV_data objects]
             base_calls_mtx: SNV-barcode matrix containing lists of base calls
             num(int): number of individual genotypes
             barcodes: list of cell barcodes
             model_genotypes(list(DataFrame): list of num model genotypes represented in snv-probdistrib DataFrame as <RR,RA,AA>
             assigned(list): lists of cell/barcode assigned to each genotype model
        """
        self.all_SNVs = all_SNVs
        self.ref_bc_mtx = base_calls_mtx[0]
        self.alt_bc_mtx = base_calls_mtx[1]
        self.barcodes = barcodes
        self.num = num
        self.model_genotypes = model_genotypes
        self.assigned = assigned

        if self.model_genotypes == []:
            for n in range(num):
                self.model_genotypes.append([])
                self.model_genotypes[n] = pd.DataFrame(
                    np.zeros((len(self.all_SNVs), 3)),
                    index=SNV_data.get_all_SNV_pos(self.all_SNVs),
                    columns=['RR', 'RA', 'AA'])

    def initialise_model_genotypes(self):
        """
	Initialises model from dirichlet distribution
	"""
        for n in range(self.num):
            dir_sim = np.random.dirichlet((1,1),len(self.all_SNVs))  # generate random dirichlet distribution to simulate genotypes probability
            gt = [item[0] for item in dir_sim]   # take the first value of each pair in the distribution results
            i = 0
            for snv in self.all_SNVs:
                pos = "{}:{}".format(snv.CHROM, snv.POS)
                self.model_genotypes[n].loc[pos] = ((1-gt[i])**2, 2*gt[i]*(1-gt[i]), gt[i]**2)   # based on Hardy-Weinberg equilibrium
                i += 1

    def calculate_model_genotypes(self):
        """
        Calculates a probability distribution <RR, RA, AA> for each SNV position based on counts of bases
        """
        coverage = [0] * self.num
        no_coverage = [0] * self.num
        for snv in self.all_SNVs:  # for each snv position
            pos = "{}:{}".format(snv.CHROM, snv.POS)
            for n in range(self.num):  # for the number of genotypes being modelled
                cluster = self.assigned[n]  # list of barcodes

                # get snv ALT base calls for cells assigned to model (cluster), returns pandas.Series
                ref_count = sum(self.ref_bc_mtx.loc[pos, cluster])
                alt_count = sum(self.alt_bc_mtx.loc[pos, cluster])
                tot_count = ref_count + alt_count
                if ref_count > 0 or alt_count > 0:  # if cells assigned in model carry reads at this snv
                    P_A_S = alt_count / tot_count       # probability of alt counts
                    coverage[n] += 1
                    self.model_genotypes[n].loc[pos] = ((1-P_A_S)**2, 2*P_A_S*(1-P_A_S), P_A_S**2)   # derive genotype probabilities using observed ALT probability, HWE

                else:
                    no_coverage[n] += 1
                    self.model_genotypes[n].loc[pos] = (BASE_RR_PROB, BASE_RA_PROB, BASE_AA_PROB)
        
        #pdb.set_trace()
        print("Coverage for model 1 = {}, no coverage = {}".format(coverage[0], no_coverage[0]))
        
        print("Coverage for model 2 = {}, no coverage = {}".format(coverage[1], no_coverage[1]))

    def assign_cells(self):
        """
        Assign cells/barcodes to most likely genotype model
        Calculate likelihood cell/barcode came from individual i (i.e. model genotype 1...n):
        For each cell/barcode c
        Product over all variants v
        Sum over all genotypes <RR,RA,AA>
        Product over all reads in cell c for variant v
        Pr(base_called | g)*P(g|s)
        where P(g|s) = set(RR,RA,AA)
        Data required:
        cell assignments : self.assigned = [[barcodes],[barcodes]]
        genotype model : self.model_genotypes
        snv pos : self.model_genotypes[0].index
        self.bc_mtx : SNV-barcode matrix containing lists of base calls
        Issue: Some barcodes in matrix have no base calls as they do not cover snp region
        Due to base_calls_mtx being instantiated from list of all barcodes in bam file (not just those covering snv)
        """
        self.old_assignment = self.assigned
        self.assigned = []
        for _ in range(self.num):
            self.assigned.append([])

        err = 0.01

        for barcode in self.barcodes:  # for each cell
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

                    ref_count = self.ref_bc_mtx.loc[pos, barcode]
                    alt_count = self.alt_bc_mtx.loc[pos, barcode]
                    tot_count = ref_count + alt_count

                    if tot_count > 0:
                        # calculate probability of alt_count in the barcode at the SNV conditioned on general genotype
                        # probability of A in RR/RA/AA is 0.01, 0.5, 0.09
                        P_A_g = binom.pmf(alt_count, tot_count, (err, 0.5, 1-err))

                    else:
                        P_A_g = (err, 0.5, 1-err)

                    # calculate P(c|S_n) for the SNV, sum of the three P_c_given_g,S
                    P_c_given_S = sum(P_A_g * self.model_genotypes[n].loc[pos])

                    # Product over all variant positions
                    all_var_llhood *= P_c_given_S

                cell_llhood.append(all_var_llhood)

            # transform into posterior probability, assuming P(S_n) are all equal (prior probability)
            cell_llhood = np.array(cell_llhood) / sum(cell_llhood)
            # if likelihoood is equal in cases such as: barcode has no coverage over any snv region
            if all(i == cell_llhood[0] for i in cell_llhood):
                pass
            else:
                n = np.argmax(cell_llhood)
                self.assigned[n].append(barcode)

        for n in range(self.num):
            self.assigned[n] = sorted(self.assigned[n])

    def compare_cell_assignments(self):
        dif = []
        for n in range(self.num):
            old = self.old_assignment[n]
            new = self.assigned[n]
            count = 0
            for item in new:
                if item in old:
                    count += 1
            if len(new) > 0:
            	print('{}% of barcodes in new also in old'.format(
                count * 100 / len(new)))
            if (len(old) + len(new)) > 0:
            	sm = jaccard_similarity(old, new)
            	dif.append(sm)
        return dif


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

def group_models(all_models):
    """Group models from multiple runs based on similarity of barcode
    assignments"""
    for e, model in enumerate(all_models):
        if e == 0:  # first iteration
            model_A = [model.assigned[0]]
            model_B = [model.assigned[1]]
        else:
            # compare current two to first model entered in A and B
            for n in range(2):
                smA = jaccard_similarity(model_A[0], model.assigned[n])
                smB = jaccard_similarity(model_B[0], model.assigned[n])
                print("model_A to model{}[{}] :".format(e + 1, n), smA)
                print("model_B to model{}[{}] :".format(e + 1, n), smB)

                match = (smA, smB).index(max(smA, smB))
                if match == 0:
                    model_A.append(model.assigned[n])
                elif match == 1:
                    model_B.append(model.assigned[n])
    return model_A, model_B


def compare_model_groups(model_group, group_name):
    for n in range(len(model_group)):
        for m in range(n + 1, len(model_group)):
            sm = jaccard_similarity(model_group[n], model_group[m])
            print("{}[{}], {}[{}]: {}".format(group_name, n + 1,
                                              group_name, m + 1, sm))


def jaccard_similarity(x, y):
    intersect = len(set.intersection(*[set(x), set(y)]))
    union = len(set.union(*[set(x), set(y)]))
    return intersect / float(union)


def run_model(all_SNVs, base_calls_mtx, barcodes, num_models):
    model = model_genotype(all_SNVs, base_calls_mtx, barcodes, num_models,
                           model_genotypes=[], assigned=None)
    dif = [0 for _ in range(num_models)]
    print("\ninitialising model genotypes", datetime.datetime.now().time())
    model.initialise_model_genotypes()
    #pdb.set_trace()
    print("initial cell assignments", datetime.datetime.now().time())
    model.assign_cells()   # initial assignment
    #pdb.set_trace()
    length = []
    for assigned in model.assigned:
        length.append(len(assigned))
    print(length)
    print("Commencing E-M")
    while any(i < 0.80 for i in
              dif):  # dif between current and prev cell assignments
        print("calculating model ", datetime.datetime.now().time())
        model.calculate_model_genotypes()
        #pdb.set_trace()
        print("assigning cells", datetime.datetime.now().time())
        model.assign_cells()
        #pdb.set_trace()
        old_dif = dif
        dif = model.compare_cell_assignments()
        length = []
        for assigned in model.assigned:
            length.append(len(assigned))
        print(length)
        print("difference:", dif)
        if dif == old_dif and sum(dif) > 0:
            break
    print("final difference:", dif)
    return model


def main():

    num_runs = 3            # number of runs
    num_models = 2          # number of models in each run


    # Mixed donor files
    bc_file = "test.txt"   # validated cell barcodes
    ref_csv = 'test_ref.csv'  # reference matrix
    alt_csv = 'test_alt.csv'  # alternative matrix
    outfile_A = 'model_A{}_de.txt'
    outfile_B = 'model_B{}_de.txt'
    
    all_SNVs = []  # list of SNV_data objects
    matrix = pd.read_csv(ref_csv)
    for record in matrix.loc[:,'SNV']:
        all_SNVs.append(SNV_data(record))
    
    print ("Starting data collection", datetime.datetime.now().time())
    
    all_POS = SNV_data.get_all_SNV_pos(all_SNVs)
    barcodes = get_all_barcodes(bc_file)
    base_calls_mtx = read_base_calls_matrix(ref_csv, alt_csv)

    # build all models of all runs (num_runs x num_models)
    all_models = []
    for run in range(num_runs):
        model = run_model(all_SNVs, base_calls_mtx, barcodes, num_models)
        all_models.append(model)
        print("Finished model {} at {}".format(run + 1,
                                               datetime.datetime.now().time()))

    print("\nFinished all models... Comparing models...")

    # separate models into groups based on similarity of barcode assignments
    model_A, model_B = group_models(all_models)

    # compare models between each grouped member
    compare_model_groups(model_A, "model_A")
    compare_model_groups(model_B, "model_B")

    # save the output barcode assignments
    for m, model in enumerate(model_A):
        file = outfile_A.format(m + 1)
        file = open(file, "w")
        for barcode in model:
            file.write(barcode + '\n')
        file.close()

    for m, model in enumerate(model_B):
        file = outfile_B.format(m + 1)
        file = open(file, "w")
        for barcode in model:
            file.write(barcode + '\n')
        file.close()


if __name__ == "__main__":
    main()

