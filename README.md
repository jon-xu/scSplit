# sc_split
### Genotype-free demultiplexing of pooled single-cell RNAseq, using a hidden state model for identifying genetically distinct samples within a mixed population.  
#### It has been used to demultiplex up to 8 samples on 10X platform.

### How to install:
    1) install python 3.6
    2) install python packages: numpy, math, pandas, pickle, pysam, random, scipy, sklearn, statistics, vcf
    3) git clone https://github.com/jon-xu/sc_split

### How to run the toolset:

##### 1. Data quality control and filtering
   *a) Copy target BAM file (barcodes marked with CB:Z: tag) into the same folder of scSplit, keep only the reads with white listed barcodes to reduce technical noises.*
   
   *b) Process BAM file from scRNA-Seq in a way that reads with any of following patterns be filtered out: quality is lower than 10,  is unmapped segment, is secondary alignment, not passing filters, is PCR or optical duplicate, or is supplementary alignment. Example: samtools view -S -b -q 10 -F 3844 original.bam > target.bam*
   
   *c) Mark BAM file for duplication, and get it sorted and indexed, using rmdup, sort, index commands in samtools*
   
##### 2. Calling for single-nucleotide variants
   *a) use freebayes v1.2 to call SNVs from the mixed sample BAM file after being processed in the first step, set the parameters for freebayes so that no insertion and deletions (indels), nor Multi-nucleotide polymorphysim (MNP) or complex events would be captured, set minimum allele count to 2 and set minimum base quality to 1.  Example: freebayes -f <reference.fa> -iXu -C 2 -q 1 target.bam snv.vcf*
   
   *b) The output VCF file will be futher filtered so that only the SNVs with quality score larger than 30 would be kept.*

##### 3. Building allele count matrices
   *a) run python script "sc_split_matrices.py" and get two .csv files ("ref_filtered.csv" and "alt_filtered.csv") as output.*

##### 4. Exectuion and verification of demultiplexing
   *a) use the two generated allele counts matrices files to demultiplex the cells into different samples.  Doublet sample will not have the same sample ID every time, which will be explicitly indicated in the log file*
   
   *b) run python script "sc_split_main.py"*
   
   *c) "sc_split_doublet.txt": indicating which cluster is doublet state*
   
   *d) "sc_split_barcodes_{n}.csv": N+1 indicating barcodes assigned to each of the N+1 samples (including doublet state)*
   
   *e) "sc_split_dist_alleles.txt": the distinguishing alleles that can be used to genotype and assign sample to clusters*
   
   *f) "sc_split_dist_matrix.csv": the ALT alelle Presence/Absence (P/A) matrix as a reference in assigning sasmple to clusters*
   
   *g) "model.found", a python pickle dump containing the final allele fraction model (model.model_MAF), and the probability of each cell belonging to each sample (model.P_s_c)*
   
   *h) "sc_split.log" log file containing information for current run, iterations, and final Maximum Likelihood and doublet sample*

##### 5. Generate genotypes based on the split result
   *a) run python script "sc_split_vcf.py"*
   
   *b) VCF file ("sc_split.vcf") will be generated for the logarithm-transformed genotype likelihoods for all sample models.*

<br/>

![alt text](https://github.com/jon-xu/sc_split/blob/master/man/figure1_pipeline.png)
