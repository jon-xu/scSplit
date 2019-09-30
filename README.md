# scSplit [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3464622.svg)](https://doi.org/10.5281/zenodo.3464622)
### Genotype-free demultiplexing of pooled single-cell RNA-seq, using a hidden state model for identifying genetically distinct samples within a mixed population.  
#### It has been tested on 3 to 8 real mixed samples using 10X pipeline, and on up to 32-mixed simulated datasets

### How to install:
  1) install python 3.6+
  2) make sure below python packages can be imported:
  
     math, numpy, pandas pickle, pysam, PyVCF, scikit-learn, scipy, statistics
  3) "git clone https://<span></span>github.com/jon-xu/scSplit" or "pip  install scSplit"
  4) run with "\<PATH\>/scSplit \<command\> \<args\>" or "python \<PATH\>/scSplit \<command\> \<args\>" 

### Overall Pipeline:

![alt text](https://github.com/jon-xu/scSplit/blob/master/man/workflow.png)

### 1. Data quality control and filtering
   a) Make sure pooled scRNA-seq BAM file doesn't contain reads from unknown barcodes, you can do this by "grep -vFwf <whitelist> <xxx>.sam > qcresult" - searching for invalid reads in SAM format of the source BAM using a file of whitelist barcodes.

   b) Filter processed BAM in a way that reads with any of following patterns be removed: read quality lower than 10,  being unmapped segment, being secondary alignment, not passing filters, being PCR or optical duplicate, or being supplementary alignment.
   
   e.g. samtools view -S -b -q 10 -F 3844 processed.bam > filtered.bam
   
   b) Mark BAM file for duplication, and get it sorted and indexed, using rmdup, sort, index commands in samtools
   
### 2. Calling for single-nucleotide variants
   a) Use freebayes v1.2 to call SNVs from the mixed sample BAM file after being processed in the first step, set the parameters for freebayes so that no insertion and deletions (indels), nor Multi-nucleotide polymorphysim (MNP) or complex events would be captured, set minimum allele count to 2 and set minimum base quality to 1.
   
   e.g. freebayes -f <reference.fa> -iXu -C 2 -q 1 filtered.bam > snv.vcf
   
   This step could take very long (up to 30 hours if not using parallel processing), GATK or other SNV calling tools should work as well.  In order to fasten the calling process, user can split the BAM by chromosome and call SNVs separately and merge the vcf files afterwards.
   
   b) The output VCF file should be futher filtered so that only the SNVs with quality score larger than 30 would be kept.
   
   c) Typical number of filtered SNVs is roughly between 20,000 and 60,000.

### 3. Building allele count matrices
   a) Run "scSplit count" and get two .csv files ("ref_filtered.csv" and "alt_filtered.csv") as output.
   
   input parameters:
      
        -v, --vcf, VCF from mixed BAM
        -i, --bam, mixed sample BAM        
        -b, --bar, barcodes whitelist
        -t, --tag, tag for barcode (default: "CB")
        -c, --com, common SNVs    
        -r, --ref, output Ref count matrix        
        -a, --alt, output Alt count matrix
        
        e.g. scSplit count -v mixed_genotype.vcf -i filtered.bam -b barcodes.tsv -r ref_filtered.csv -a alt_filtered.csv
   
   b) It is **strongly recommended** to use below SNV list to filter the matrices to improve prediction accuracy:

      Common SNPs (e.g. Human common SNPs from 1000 Genome project)
   
      hg19: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
   
      hg38: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/
        
      To process the genotype files of common SNPs, either download per-chromosome files and concatenate them using bcftools or download the whole genome file, take the first two columns of the vcf file and replace the tab with colon sign so that each line is one SNV, e.g., "1:10177". 
      
      Processed common SNVs for hg19 and hg38 can be found here: http://data.genomicsresearch.org/Projects/scSplit/CommonSNVs

   Please specify the common SNVs in scSplit count using -c/--com parameter, please make sure your common SNVs list does not have header row.
   
   c) This step could be memory consuming, if the number of SNVs and/or cells are high. As a guideline, building matrices for 60,000 SNVs and 10,000 cells might need more than 30GB RAM to run, please allow enough RAM resource for running the script.
   
   d) Typical runtime for this step is about one hour, depending on the nature of the data and the resources being allocated.

### 4. Demultiplexing and generate ALT P/A matrix
   a) Use the two generated allele counts matrices files to demultiplex the cells into different samples.  Doublet sample will not have the same sample ID every time, which will be explicitly indicated in the log file

   b) Run "scSplit run" with input parameters:
      
        -r, --ref, input Ref count matrix        
        -a, --alt, input Alt count matrix        
        -n, --num, expected number of mixed samples (-n 0: autodetect mode)
        -s, --sub, (optional) maximum number of subpopulations in autodetect mode, default: 10
        -e, --ems, (optional) number of EM repeats to avoid local maximum, default: 30
        -d, --dbl, (optional) correction for doublets, 0 for no doublets, and no refinement on the results if not specified or specified percentage is less than detected
        -v, --vcf, (optional) known individual genotypes to map clusters and samples using distinguishing variants

        e.g. scSplit run -r ref_filtered.csv -a alt_filtered.csv -n 8
        
        # below command will tell the script to expect 20% doublets if the natually found doublets are less than that:
        e.g. scSplit run -r ref_filtered.csv -a alt_filtered.csv -n 8 -d 0.2
        
        # (beta) -n 0 -s <sub>, let system decide the optimal sample number between 2 and <sub>
        e.g. scSplit run -r ref_filtered.csv -a alt_filtered.csv -n 0 -s 12

   c) Below files will be generated:

      "scSplit_result.csv": barcodes assigned to each of the N+1 cluster (N singlets and 1 doublet cluster), doublet marked as DBL-<n> (n stands for the cluster number)
      "scSplit_dist_variants.txt": the distinguishing variants that can be used to genotype and assign sample to clusters
      "scSplit_dist_matrix.csv": the ALT allele Presence/Absence (P/A) matrix on distinguishing variants for all samples as a reference in assigning sample to clusters, NOT including the doublet cluster, whose sequence number would be different every run (please pay enough attention to this)
      "scSplit_PA_matrix.csv": the full ALT allele Presence/Absence (P/A) matrix for all samples, NOT including the doublet cluster, whose sequence number would be different every run (please pay enough attention to this)
      "scSplit_P_s_c.csv", the probability of each cell belonging to each sample
      "scSplit.log" log file containing information for current run, iterations, and final Maximum Likelihood and doublet sample
      
   d) This step is also memory consuming, and the RAM needed is highly dependent on the quantity of SNVs from last step and the number of cells. As a guideline, a matrix with 60,000 SNVs and 10,000 cells might need more than 50GB RAM to run, please allow enough RAM resource for running the script.
   
   e) Typical runtime for this step is about half an hour, with default parameters, depending on the nature of the data and the resources being allocated.

### 5. (Optional) Generate sample genotypes based on the split result
   a) Run "scSplit genotype" with input parameters:
       
        -r, --ref, Ref count CSV as output        
        -a, --alt, Alt count CSV as output
        -p, --psc, generated P(S|C)

        e.g. scSplit genotype -r ref_filtered.csv -a alt_filtered.csv -p scSplit_P_s_c.csv
        
   b) VCF file ("scSplit.vcf") will be generated for the logarithm-transformed genotype likelihoods for all sample models.

<br/>

