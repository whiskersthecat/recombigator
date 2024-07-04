# recombigator
a tool to classify Nanopore or HiFi DNA reads as belonging to genome A, genome B, or recombinant

## Installation
**recombigator** is compiled on the command line using:
> g++ recombigator.cc -std=c++11 -O3 -o recombigator

## Usage
> ./recombigator alignment_A.sam alignment_B.sam variations.tsv [output_name]

## Input Files
### alignment_A.sam
includes all reads aligned to genome A, with match/mismatch information in the CIGAR string

for example, to align using the minimap2 aligner (https://github.com/lh3/minimap2):
> ./minimap2 --eqx -a --secondary=no genomeA.fasta reads.fasta

### variations.tsv
a tab seperated value file with the format specified by SYRI (https://schneebergerlab.github.io/syri/fileformat.html):
> genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ 
**region_type** beginning with "SYN" and **variation_type** = "SNP" indicates a SNP in a syntenic region

**genomeA_chrom** and **genomeA_pos** give the coordinates of the SNP in genome A

fields marked with ~ are not used and can be any string

## Output
recombigator classifies by reads by looking at the SNPs spanned by the alignments

in order to partially mitigate bias caused by sequencing errors, SNPs with mismatches in both alignments are ignored

each read is sorted into one of the following categories:
- **GENOTYPEA**
- **GENOTYPEB**
- **PERIPHERAL_RECOMBINANT** recombinant with low confidence (only one SNP for one genotype)
- **RECOMBINANT**
- **CHOMP** reads can be chomped because:
1. Aligned to different chromosomes
2. Consensus SNP Matrix switches gt more than once
3. Less than two SNPs in Consensus SNP Matrix
