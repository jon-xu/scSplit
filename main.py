"""
Reference free AF-based demultiplexing on pooled scRNA-seq (state intialisation using pca, with random factor to avoid local maximum)
Jon Xu (jun.xu@uq.edu.au)
Lachlan Coin
Aug 2018
"""

import numpy as np
import pandas as pd
import vcf
import statistics as stat
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import datetime, pickle, argparse

class models:
    def __init__(self, base_calls_mtx, num):
        """
        Model class containing SNVs, matrices counts, barcodes, model allele fraction with assigned cells

        Parameters:
             base_calls_mtx(list of Dataframes): SNV-barcode matrix containing lists of base calls
             ref_bc_mtx: SNV-barcode matrix for reference allele counts
             alt_bc_mtx: SNV-barcode matrix for alternative allele counts
             all_POS(list): list of SNVs positions
             barcodes(list): list of cell barcodes
             num(int): number of total samples
             P_s_c(DataFrame): barcode/sample matrix containing the probability of seeing sample s with observation of barcode c
             lP_c_s(DataFrame): barcode/sample matrix containing the log likelihood of seeing barcode c under sample s, whose sum should increase by each iteration
             assigned(list): final lists of cell/barcode assigned to each cluster/model
             model_af(Dataframe): Dataframe of model allele frequencies P(A) for each SNV and state
        """
        self.ref_bc_mtx, self.alt_bc_mtx = base_calls_mtx[0], base_calls_mtx[1]
        self.all_POS, self.barcodes = base_calls_mtx[2].tolist(), base_calls_mtx[3].tolist()
        self.num = int(num) + 1  # additional doublet state
        self.P_s_c = pd.DataFrame(0, index = self.barcodes, columns = range(self.num))
        self.lP_c_s = pd.DataFrame(0, index = self.barcodes, columns = range(self.num))
        self.assigned, self.reassigned = [], []
        for _ in range(self.num):
            self.assigned.append([])
            self.reassigned.append([])
        self.model_af = pd.DataFrame(0, index=self.all_POS, columns=range(self.num))
        self.pseudo = 1

        # set background alt count proportion as allele fraction for each SNVs of doublet state, with pseudo count added for 0 counts on multi-base SNPs
        N_A = self.alt_bc_mtx.sum(axis=1) + self.pseudo
        N_R = self.ref_bc_mtx.sum(axis=1) + self.pseudo
        N_T = N_A + N_R
        k_ref,  k_alt = N_R / N_T, N_A / N_T

        # find barcodes for state initialization, using subsetting/PCA/K-mean
        base_mtx = (self.alt_bc_mtx + self.ref_bc_mtx).toarray()
        rows, cols = [*range(base_mtx.shape[0])], [*range(base_mtx.shape[1])]
        irows, icols = np.array(rows), np.array(cols)
        nrows, ncols = len(rows), len(cols)
        mrows = mcols = 0
        while (mrows < 0.9 * nrows or mcols < 0.9 * ncols) and (nrows >= 10 or ncols >= 10):
            rbrows = np.sort(np.unique(list(map(int, np.random.beta(1,10,int(0.1*nrows))*nrows))))
            rbcols = np.sort(np.unique(list(map(int, np.random.beta(1,10,int(0.1*ncols))*ncols))))
            rows = np.count_nonzero(base_mtx, axis=1).argsort().tolist()
            cols = np.count_nonzero(base_mtx, axis=0).argsort().tolist()
            rows = [item for index, item in enumerate(rows) if index not in rbrows]
            cols = [item for index, item in enumerate(cols) if index not in rbcols]
            irows, icols = irows[rows], icols[cols]
            nrows, ncols = len(rows), len(cols)
            base_mtx = base_mtx[rows][:,cols]
            mrows = min(np.count_nonzero(base_mtx, axis=0))
            mcols = min(np.count_nonzero(base_mtx, axis=1))
        alt_subset = self.alt_bc_mtx[irows][:, icols].todense()
        ref_subset = self.ref_bc_mtx[irows][:, icols].todense()
        alt_prop = (alt_subset + 0.01) / (alt_subset + ref_subset + 0.02)
        alt_pca = StandardScaler().fit_transform(alt_prop.T)
        pca = PCA(n_components=min(nrows, ncols, 20))
        pca_alt = pca.fit_transform(alt_pca)
        if pca_alt.shape[0] < self.num:
            print('not enough informative cells to support model initialization')
        else:
            kmeans = KMeans(n_clusters=self.num, random_state=0).fit(pca_alt)
            # intialise allele frequency for model states
            self.initial = []
            for n in range(self.num):
                self.initial.append([])
                for index, col in enumerate(icols):
                    if kmeans.labels_[index] == n:
                        self.initial[n].append(self.barcodes[col])
                barcode_alt = np.array(self.alt_bc_mtx[:, icols[kmeans.labels_==n]].sum(axis=1))
                barcode_ref = np.array(self.ref_bc_mtx[:, icols[kmeans.labels_==n]].sum(axis=1))
                self.model_af.loc[:, n] = (barcode_alt + k_alt) / (barcode_alt + barcode_ref + k_alt + k_ref)


    def run_EM(self):

        # commencing E-M
        iterations = 0
        self.sum_log_likelihood = [1,2]
        while self.sum_log_likelihood[-2] != self.sum_log_likelihood[-1]:
            iterations += 1
            progress = 'Iteration ' + str(iterations) + '   ' + str(datetime.datetime.now()) + '\n'
            with open('scSplit.log', 'a') as myfile: myfile.write(progress)
            self.calculate_cell_likelihood()  # E-step
            self.calculate_model_af()  # M-step
            self.sum_log_likelihood.append(self.lP_c_m)


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
            self.lP_c_s.loc[:, n] = matcalc.sum(axis=0).tolist()[0]

        # transform to cell sample probability P(s|c) using Baysian rule
        for i in range(self.num):
            denom = 0
            for j in range(self.num):
                denom += 2 ** (self.lP_c_s.loc[:, j] - self.lP_c_s.loc[:, i])
            self.P_s_c.loc[:, i] = 1 / denom

        # calculate model likelihood: logP(x|theta) = log{Sigma_yP(x,y|theta)}
        self.lP_c_m = ((self.lP_c_s.subtract(self.lP_c_s.min(axis=1), axis=0).pow(2).sum(axis=1) + 1).apply(np.log2) + self.lP_c_s.min(axis=1)).sum(axis=0)


    def calculate_model_af(self):
        """
        Update the model allele fraction by distributing the alt and total counts of each barcode on a certain snv to the model based on P(s|c)
        """

        N_ref = self.ref_bc_mtx.sum(axis=1) + self.pseudo
        N_alt = self.alt_bc_mtx.sum(axis=1) + self.pseudo
        k_ref, k_alt = N_ref / (N_ref + N_alt), N_alt / (N_ref + N_alt)
        self.model_af = pd.DataFrame((self.alt_bc_mtx.dot(self.P_s_c) + k_alt) / ((self.alt_bc_mtx + self.ref_bc_mtx).dot(self.P_s_c) + k_ref + k_alt),
                                        index = self.all_POS, columns = range(self.num))


    def assign_cells(self):
        """
            Final assignment of cells according to P(s|c) >= 0.9
        """

        for n in range(self.num):
            self.assigned[n] = sorted(self.P_s_c.loc[self.P_s_c[n] >= 0.99].index.values.tolist())


    def define_doublet(self):
        """
            Locate the doublet state
        """
        cross_state = pd.DataFrame(0, index = range(self.num), columns = range(self.num))
        for i in range(self.num):
            for j in range(self.num):
                index = []
                # transform barcode assignments to indices
                for item in self.assigned[j]:
                    index.append(self.barcodes.index(item))
                matcalc = self.alt_bc_mtx.T.multiply(self.model_af.loc[:, i].apply(np.log2)).T \
                        + self.ref_bc_mtx.T.multiply((1 - self.model_af.loc[:, i]).apply(np.log2)).T
                cross_state.loc[j, i] = matcalc[:, index].sum()
        result = cross_state.sum(axis=0).tolist()
        self.doublet = result.index(max(result))


    def refine_doublets(self):
        """
            Find falsely assigned doublets
        """
        N_ref_mtx, N_alt_mtx = pd.DataFrame(0, index=self.all_POS, columns=range(self.num)), pd.DataFrame(0, index=self.all_POS, columns=range(self.num))
        found = []
        self.reassigned[self.doublet] = self.assigned[self.doublet]
        for n in range(self.num):
            bc_idx = [i for i, e in enumerate(self.barcodes) if e in self.assigned[n]]
            # REF/ALT alleles counts from cells assigned to state n
            N_ref_mtx.loc[:, n], N_alt_mtx.loc[:, n] = self.ref_bc_mtx[:, bc_idx].sum(axis=1), self.alt_bc_mtx[:, bc_idx].sum(axis=1)
        read_depth = self.alt_bc_mtx.T.dot(np.int64((N_alt_mtx > 0) | (N_ref_mtx > 0))) * (self.P_s_c >= 0.99)
        for n in range(self.num):
            if n != self.doublet:
                count = np.count_nonzero(read_depth[n])
                found = read_depth.index[read_depth[n].argsort()[int((len(read_depth.index) - count * 0.05)) : len(read_depth.index)]]
                self.reassigned[self.doublet] += found.tolist()
                self.reassigned[n] = [x for x in self.assigned[n] if x not in found]


    def distinguishing_alleles(self, pos=[]):
        """
            Locate the distinguishing alleles
            N_ref_mtx, N_alt_mtx: SNV-state matrix for ref/alt counts in each state
        """
        
        # build SNV-state matrices for ref and alt counts
        self.dist_variants, todo = [], []
        if len(pos) != 0:
            snv = [self.all_POS[i] for i in pos]
            N_ref_mtx, N_alt_mtx = pd.DataFrame(0, index=snv, columns=range(self.num)), pd.DataFrame(0, index=snv, columns=range(self.num))
        else:
            N_ref_mtx, N_alt_mtx = pd.DataFrame(0, index=self.all_POS, columns=range(self.num)), pd.DataFrame(0, index=self.all_POS, columns=range(self.num))

        for n in range(self.num):
            bc_idx = [i for i, e in enumerate(self.barcodes) if e in self.reassigned[n]]
            # REF/ALT alleles counts from cells assigned to state n
            if len(pos) == 0:
                N_ref_mtx.loc[:, n], N_alt_mtx.loc[:, n] = self.ref_bc_mtx[:, bc_idx].sum(axis=1), self.alt_bc_mtx[:, bc_idx].sum(axis=1)
            else:
                N_ref_mtx.loc[:, n], N_alt_mtx.loc[:, n] = self.ref_bc_mtx[pos][:, bc_idx].sum(axis=1), self.alt_bc_mtx[pos][:, bc_idx].sum(axis=1)

        # judge N(A) or N(R) for each cluster
        alt_or_ref = ((N_alt_mtx >= 10) * 1 - ((N_ref_mtx >= 10) & (N_alt_mtx == 0)) * 1).drop(self.doublet, axis=1).astype(np.int8)
        alt_or_ref[alt_or_ref == 0], alt_or_ref[alt_or_ref == -1] = float('NaN'), 0     # formatting data for further analysis
        alt_or_ref = alt_or_ref.ix[[x for x in alt_or_ref.index if x[0] not in ['X','Y','MT']]]

        start, ncols = 0, self.num - 1
        # find unique alleles for each column window and then merge
        while start < ncols - 1:                          
            least_ones = []
            while len(least_ones) < (ncols - start):
                selected = []                        
                submatrix = alt_or_ref.iloc[:, start:ncols]
                # informative alleles for the selected clusters with no NAs
                if start < ncols: informative_sub = submatrix[(submatrix.var(axis=1) > 0) & (submatrix.count(axis=1) == submatrix.shape[1])]
                else: informative_sub = submatrix[submatrix.count(axis=1) == submatrix.shape[1]]
                if informative_sub.index.values.size > 0:
                    patt = informative_sub.astype(str).values.sum(axis=1)
                    unq = np.unique(patt, return_inverse=True)                  
                    for j in range(len(unq[0])):
                        # first informative alleles for each unique pattern in the cluster screen which has maximum non-NA values in original information matrix
                        selected.append(alt_or_ref.loc[informative_sub.iloc[[i for i, x in enumerate(unq[1]) if x == j]].index].count(axis=1).idxmax())
                    subt = alt_or_ref.loc[np.unique(selected)].iloc[:, start:ncols] 
                    d = np.linalg.svd(subt, full_matrices=False)[1]
                    if sum(d > 1e-10) >= (ncols - start):
                        least_ones = [subt[subt.sum(axis=1) == min(subt.sum(axis=1))].index[0]]
                        while len(least_ones) < (ncols - start):            # Gram-Schmidt Process
                            svd = np.linalg.svd(subt.loc[least_ones].transpose(), full_matrices=False)
                            U, V = svd[0], svd[2]
                            if len(least_ones) == 1: D=np.asmatrix(svd[1])
                            else: D = np.diag(svd[1])
                            colU, Vinv, Dinv = U.shape[1], np.linalg.solve(V, np.diag([1]*len(V))), np.linalg.solve(D, np.diag([1]*len(D)))
                            proj = np.asmatrix(np.zeros((colU,subt.shape[0])))
                            for j in range(subt.shape[0]):
                                for i in range(colU):
                                    proj[i, j] = np.matmul(U[:,i], subt.iloc[j])
                            W = np.matmul(np.matmul(Vinv, Dinv), proj)
                            R = np.matmul(subt.loc[least_ones].transpose(), W)
                            diff = (subt.transpose() - R).transpose()
                            subt1 = subt.loc[diff[diff.var(axis=1) > (0.5 * max(diff.var(axis=1)))].index]
                            least_ones += [subt1[subt1.sum(axis=1) == min(subt1.sum(axis=1))].index[0]]
                ncols -= 1
            self.dist_variants += least_ones
            start = ncols + 1
            ncols = self.num - 1

        self.dist_variants = list(set(self.dist_variants))

        # number of rows to distinguish each cluster pair
        col_diff = pd.DataFrame(0, index=alt_or_ref.columns, columns=alt_or_ref.columns)
        for i in range(1, col_diff.shape[0]):
            for j in range(i):
                # number of alleles that cluster i and j are different in
                col_diff.iloc[i, j] = (alt_or_ref.loc[self.dist_variants].iloc[:, [i, j]].var(axis=1) > 0).sum()
                if col_diff.iloc[i, j] == 0:
                    todo.append([i, j])

        # for non-covered pairs, expand to get distinguishing alleles
        for pair in todo:
            try:
                new =alt_or_ref[alt_or_ref.loc[:, pair].var(axis=1) > 0].index[0]
                self.dist_variants.append(new)
            except:
                with open('scSplit.log', 'a') as myfile: myfile.write('\n not all distinguish-able! \n')

        self.dist_variants = list(set(self.dist_variants))
        self.dist_matrix = alt_or_ref.loc[self.dist_variants]


def main():

    # Process command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', required=True,  help='Ref count CSV input')
    parser.add_argument('-a', '--alt', required=True,  help='Alt count CSV input')
    parser.add_argument('-n', '--num', required=True,  help='Number of mixed samples')
    parser.add_argument('-v', '--vcf', required=False, help='vcf file for filtering distiinguishing variants')
    args = parser.parse_args()

    progress = 'Starting data collection: ' + str(datetime.datetime.now()) + '\n'
    with open('scSplit.log', 'a') as myfile: myfile.write(progress)

    # Read in existing matrix from the csv files
    ref = pd.read_csv(args.ref, header=0, index_col=0)
    alt = pd.read_csv(args.alt, header=0, index_col=0)
    ref_s, alt_s = csr_matrix(ref.values), csr_matrix(alt.values)
    base_calls_mtx = [ref_s, alt_s, ref.index, ref.columns]
    progress = 'Allele counts matrices uploaded: ' + str(datetime.datetime.now()) + '\n'
    with open('scSplit.log', 'a') as myfile: myfile.write(progress)

    max_likelihood = -1e10
    for i in range(30):
        with open('scSplit.log', 'a') as myfile: myfile.write('round ' + str(i) + '\n')

        model = models(base_calls_mtx, args.num)
        if model.model_af.sum().sum() > 0:
            model.run_EM()
            model.assign_cells()
            if model.sum_log_likelihood[-1] > max_likelihood:
                max_likelihood = model.sum_log_likelihood[-1]
                initial, assigned, af, p_s_c = model.initial, model.assigned, model.model_af, model.P_s_c
    model.assigned, model.initial, model.model_af, model.P_s_c = assigned, initial, af, p_s_c
    model.define_doublet()
    model.refine_doublets()

    pos = []
    try:
        for record in vcf.Reader(open(args.vcf, 'r')):
            # only keep high R2 variants
            try:
                if float(record.INFO['R2'][0]) > 0.9:
                    pos.append(model.all_POS.index(str(record.CHROM)+':'+str(record.POS)))
            except:
                continue
    except:
        pass
        
    model.distinguishing_alleles(pos)

    # generate outputs
    with open('scSplit_model', 'wb') as f:
        pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)
    for n in range(int(args.num)+1):
        with open('scSplit_barcodes_{}.csv'.format(n), 'w') as myfile:
            for item in model.reassigned[n]:
                myfile.write(str(item) + '\n')
    model.P_s_c.to_csv('scSplit_P_s_c.csv')
    with open('scSplit_dist_variants.txt', 'w') as myfile:
        for item in model.dist_variants:
            myfile.write(str(item) + '\n')
    with open('scSplit_doublet.txt', 'w') as myfile:
        myfile.write('Cluster ' + str(model.doublet) + ' is doublet.\n')
    model.dist_matrix.to_csv('scSplit_dist_matrix.csv')


if __name__ == '__main__':
    main()
