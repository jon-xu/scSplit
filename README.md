# scSplit
### Genotype-free demultiplexing of pooled single-cell RNA-seq, using a hidden state model for identifying genetically distinct samples within a mixed population.  
#### It has been used to demultiplex 3 to 8 samples on 10X datasets.

### How to install:
  1) install python 3.6+
  2) make sure below python packages can be imported:
  
     numpy, math, pandas, pickle, pysam, random, scipy, statistics, scikit-learn, PyVCF
  3) git clone https://github.com/jon-xu/scSplit or pip install scSplit

### How to run the toolset:

![alt text](https://github.com/jon-xu/scSplit/blob/master/man/workflow.png)

### 1. Data quality control and filtering
   a) Copy target BAM file (barcodes marked with CB:Z: tag) into the same folder of scSplit, keep only the reads with white listed barcodes to reduce technical noises.
   
   b) Process BAM file from scRNA-Seq in a way that reads with any of following patterns be filtered out: quality is lower than 10,  is unmapped segment, is secondary alignment, not passing filters, is PCR or optical duplicate, or is supplementary alignment.
   
   E.g.: samtools view -S -b -q 10 -F 3844 original.bam > target.bam
   
   c) Mark BAM file for duplication, and get it sorted and indexed, using rmdup, sort, index commands in samtools
   
### 2. Calling for single-nucleotide variants
   a) Use freebayes v1.2 to call SNVs from the mixed sample BAM file after being processed in the first step, set the parameters for freebayes so that no insertion and deletions (indels), nor Multi-nucleotide polymorphysim (MNP) or complex events would be captured, set minimum allele count to 2 and set minimum base quality to 1.
   
   E.g.: freebayes -f <reference.fa> -iXu -C 2 -q 1 target.bam > snv.vcf
   
   This step could take very long (up to 30 hours if not using parallel processing), GATK or other SNV calling tools might work as well.  Users can also split the BAM by chromosome and call SNVs separately.
   
   b) The output VCF file will be futher filtered so that only the SNVs with quality score larger than 30 would be kept.

### 3. Building allele count matrices
   a) Run python script "matrices.py" and get two .csv files ("ref_filtered.csv" and "alt_filtered.csv") as output.
      input parameters:
      
        -v, --vcf, VCF from mixed BAM
        -i, --bam, mixed sample BAM        
        -b, --barcodes, barcodes whitelist        
        -r, --ref, Ref count CSV as output        
        -a, --alt, Alt count CSV as output
        
   E.g.: python matrices.py -v mixed_genotype.vcf -i mixed.bam -b barcodes.tsv -r ref_filtered.csv -a alt_filtered.csv
   
   b) This step is memory consuming, and the RAM needed is highly dependent on the quantity of SNVs from last step and the number of cells. As a guideline, a matrix with 60,000 SNVs and 10,000 cells might need more than 30GB RAM to run, please allow enough RAM resource for running the script.

   c) Two filtering methods can be used to filter the matrices to improve prediction accuracy:
      1) Common SNPs (e.g., Human common SNPs from 1000 Genome project)
   
      hg19: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
   
      hg38: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/
        
      To process the genotype files of common SNPs, either download per-chromosome files and concatenate them using bcftools or the whole genome file, take the first two columns of the vcf file and replace the tab with colon sign so that each line is one SNV, e.g., "1:10177". Processed common SNVs for hg19 and hg38 can be found here: http://data.genomicsresearch.org/Projects/scSplit/

      2) repeat sequence regions (repeatmasker)

      hg19: http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz

      hg38: http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz

      To process the files, take the chromosome, start and end columns of the files, and filter the called SNVs with these blacklist regions and get a list of SNVs without falling into these regions.

   Then filter the matrices generated in the last step ("ref_filtered.csv" and "alt_filtered.csv") with the list of common SNVs and use them as reference and alternative matrices as inputs for scSplit run.


### 4. Demultiplexing and generate ALT P/A matrix
   a) Use the two generated allele counts matrices files to demultiplex the cells into different samples.  Doublet sample will not have the same sample ID every time, which will be explicitly indicated in the log file

   b) This step is also memory consuming, and the RAM needed is highly dependent on the quantity of SNVs from last step and the number of cells. As a guideline, a matrix with 60,000 SNVs and 10,000 cells might need more than 50GB RAM to run, please allow enough RAM resource for running the script.
   
   c) Run python script "main.py"
      input parameters:
      
        -r, --ref, Ref count CSV as input        
        -a, --alt, Alt count CSV as input        
        -n, --num, Number of mixed samples
        -v, --vcf, individual genotypes to check distinguishing variants against (optional)
        
   E.g.: python main.py -r ref_filtered.csv -a alt_filtered.csv -n 8
        
   d) "scSplit_doublet.txt": indicating which cluster is doublet state
   
   e) "scSplit_barcodes_{n}.csv": N+1 indicating barcodes assigned to each of the N+1 samples (including doublet state)
   
   f) "scSplit_dist_alleles.txt": the distinguishing alleles that can be used to genotype and assign sample to clusters
   
   g) "scSplit_dist_matrix.csv": the ALT allele Presence/Absence (P/A) matrix on distinguishing variants for all samples as a reference in assigning sample to clusters

   h) "scSplit_PA_matrix.csv": the full ALT allele Presence/Absence (P/A) matrix for all samples
   
   i) "scSplit.model", a python pickle dump containing the final allele fraction model (model.model_MAF), and the probability of each cell belonging to each sample (model.P_s_c)
   
   j) "scSplit.log" log file containing information for current run, iterations, and final Maximum Likelihood and doublet sample

### 5. Generate sample genotypes based on the split result
   a) Run python script "genotype.py"
       input parameters:
       
        -r, --ref, Ref count CSV as output        
        -a, --alt, Alt count CSV as output
        -p, --psc, generated P(S|C)
        
   E.g.: python genotype.py -r ref_filtered.csv -a alt_filtered.csv -p scSplit_P_s_c.csv
        
   b) VCF file ("scSplit.vcf") will be generated for the logarithm-transformed genotype likelihoods for all sample models.

<br/>

