# scSplit
### Genotype-free demultiplexing of pooled single-cell RNAseq, using a hidden state model for identifying genetically distinct samples within a mixed population.  
#### It has been used to demultiplex up to 8 samples on 10X platform.

### How to install:
    1) install python 3.6+
    2) install python packages: numpy, math, pandas, pickle, pysam, random, scipy, scikit-learn, PyVCF
    3) pip install scSplit -> python3 -m scSplit ..., or: git clone https://github.com/jon-xu/scSplit

### How to run the toolset:

##### 1. Data quality control and filtering
   *a) Copy target BAM file (barcodes marked with CB:Z: tag) into the same folder of scSplit, keep only the reads with white listed barcodes to reduce technical noises.*
   
   *b) Process BAM file from scRNA-Seq in a way that reads with any of following patterns be filtered out: quality is lower than 10,  is unmapped segment, is secondary alignment, not passing filters, is PCR or optical duplicate, or is supplementary alignment. Example: samtools view -S -b -q 10 -F 3844 original.bam > target.bam*
   
   *c) Mark BAM file for duplication, and get it sorted and indexed, using rmdup, sort, index commands in samtools*
   
##### 2. Calling for single-nucleotide variants
   *a) Use freebayes v1.2 to call SNVs from the mixed sample BAM file after being processed in the first step, set the parameters for freebayes so that no insertion and deletions (indels), nor Multi-nucleotide polymorphysim (MNP) or complex events would be captured, set minimum allele count to 2 and set minimum base quality to 1.  Example: freebayes -f <reference.fa> -iXu -C 2 -q 1 target.bam snv.vcf. This step could take very long (up to 30 hours if not using parallel processing), GATK or other SNV calling tools might work as well.  Users can also split the BAM by chromosome and call SNVs separately.*
   
   *b) The output VCF file will be futher filtered so that only the SNVs with quality score larger than 30 would be kept.*

##### 3. Building allele count matrices
   *a) Run python script "matrices.py" and get two .csv files ("ref_filtered.csv" and "alt_filtered.csv") as output.*
   
        -v, --vcf, VCF from mixed BAM
        -i, --bam, mixed sample BAM        
        -b, --barcodes, barcodes whitelist        
        -r, --ref, Ref count CSV as output        
        -a, --alt, Alt count CSV as output
   
   *b) This step is memory consuming, and the RAM needed is highly dependent on the quantity of SNVs from last step and the number of cells. As a guideline, a matrix with 60,000 SNVs and 10,000 cells might need more than 30GB RAM to run, please allow enough RAM resource for running the script.

##### 4. Exectuion and verification of demultiplexing
   *a) Use the two generated allele counts matrices files to demultiplex the cells into different samples.  Doublet sample will not have the same sample ID every time, which will be explicitly indicated in the log file*

   *b) This step is also memory consuming, and the RAM needed is highly dependent on the quantity of SNVs from last step and the number of cells. As a guideline, a matrix with 60,000 SNVs and 10,000 cells might need more than 50GB RAM to run, please allow enough RAM resource for running the script.
   
   *b) Run python script "main.py"*
   
        -r, --ref, Ref count CSV as input        
        -a, --alt, Alt count CSV as input        
        -n, --num, Number of mixed samples
        -v, --vcf, individual genotypes to check distinguishing variants against (optional)
        
   *c) "scSplit_doublet.txt": indicating which cluster is doublet state*
   
   *d) "scSplit_barcodes_{n}.csv": N+1 indicating barcodes assigned to each of the N+1 samples (including doublet state)*
   
   *e) "scSplit_dist_alleles.txt": the distinguishing alleles that can be used to genotype and assign sample to clusters*
   
   *f) "scSplit_dist_matrix.csv": the ALT allele Presence/Absence (P/A) matrix on distinguishing variants for all samples as a reference in assigning sample to clusters*

   *g) "scSplit_PA_matrix.csv": the full ALT allele Presence/Absence (P/A) matrix for all samples* 
   
   *h) "scSplit.model", a python pickle dump containing the final allele fraction model (model.model_MAF), and the probability of each cell belonging to each sample (model.P_s_c)*
   
   *i) "scSplit.log" log file containing information for current run, iterations, and final Maximum Likelihood and doublet sample*

##### 5. Generate genotypes based on the split result
   *a) Run python script "genotype.py"*
             
        -r, --ref, Ref count CSV as output        
        -a, --alt, Alt count CSV as output
        -p, --psc, generated P(S|C)
        
   *b) VCF file ("scSplit.vcf") will be generated for the logarithm-transformed genotype likelihoods for all sample models.*

<br/>

![alt text](https://github.com/jon-xu/scSplit/blob/master/man/figure1_pipeline.png)
