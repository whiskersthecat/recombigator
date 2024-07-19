# recombigator
A tool to classify Nanopore or HiFi DNA reads as belonging to genome A, genome B, or recombinant.

## Installation
**recombigator** is compiled on the command line using:
> g++ recombigator.cc -std=c++11 -O3 -o recombigator

## Usage
> ./recombigator alignment_A.sam alignment_B.sam syri.out [output_name]

## Input Files
### alignment_A.sam
All reads aligned to genome A, with match/mismatch information in the CIGAR string. 
For example, to align using the minimap2 aligner (https://github.com/lh3/minimap2):
> ./minimap2 --eqx -a --secondary=no genomeA.fasta reads.fasta > alignment_A.sam

### syri.out
Tab seperated value file with the format specified by SYRI (https://schneebergerlab.github.io/syri/fileformat.html):
> genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ 

Where **region_type** beginning with "SYN" and **variation_type** = "SNP" indicates a SNP in a syntenic region. 
**genomeA_chrom** and **genomeA_pos** give the coordinates of the SNP in genome A. 
Fields marked with ~ are not used and can be any string.

## Output
Recombigator classifies by reads by looking at the SNPs spanned by the alignments.
Each read is sorted into one of the following categories:
- **GENOTYPE A**
- **GENOTYPE B**
- **LOW QUALITY RECOMBINANT** recombinant with only one SNP marked as one genotype
- **MEDIUM QUALITY RECOMBINANT** recombinant with only two SNPs marked as one genotype
- **HIGH QUALITY RECOMBINANT** other recombinants
- **CHOMP : BAD ALIGN** Aligned to different chromosomes
- **CHOMP : BAD GT** Consensus SNP Matrix switches gt more than once
- **CHOMP : NO SNPS** Less than two SNPs in Consensus SNP Matrix (e.g. read not from syntenic region)
