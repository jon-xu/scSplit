# sc_split
### Genotype-free demultiplexing of pooled single-cell RNAseq, using a hidden state model for identifying genetically distinct samples within a mixed population.  It has been used to demultiplex up to 8 samples on 10X platform.

### General pipeline:

##### 1. Data quality control and filtering
   1) process BAM file from scRNA-Seq in a way that reads with any of following patterns be filtered out: quality is lower than 10,  is unmapped segment, is secondary alignment, not passing filters, is PCR or optical duplicate, or is supplementary alignment.  The BAM file will be further filtered to keep only the reads with white listed barcodes to reduce technical noises.  Finally it will be marked for duplication, sorted and indexed.
2. Calling for single-nucleotide variants, use freebayes v1.2 \cite{freebayes} to call SNVs from the mixed sample BAM file after being processed in the first step, set the parameters for freebayes so that no insertion and deletions (indels), nor Multi-nucleotide polymorphysim (MNP),  or complex events would be captured, set minimum allele count to 2 and set minimum base quality to 1.  The output VCF file will be futher filtered so that only the SNVs with quality score larger than 30 would be kept. 
3. Building allele count matrices, run the "build\_matrices.py" script and get two .csv files as output. 
4. Exectuion and verification of demultiplexing, use the two generated allele counts matrices files to demultiplex the cells into their original sample models.  There will be n+3 files for n samples as output: n .csv files to list all the barcodes assigned to all n samples, one file to capture doublet barcodes, one file to record the final allele fraction model, and the last file to record the probability of each cell linked with all sample models. 
5. Output of Genotype likelihoods for models, a VCF file will be generated for the logarithm-transformed genotype likelihoods for all sample models. 

### How to install:

### How to Run:

### How to interpret the results:
