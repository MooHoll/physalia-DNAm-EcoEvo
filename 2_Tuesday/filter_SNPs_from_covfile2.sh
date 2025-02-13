#!/bin/bash

# Use bedtools to filter a CpG-merged bismark cov file directly for SNPs.
# This accounts for SNPs at the neighbouring C on the opposite strand.
# requires bedtools (should be installed in conda environment already)

# usage: sh filter_SNPs_from_covfile2.sh <input.cov.gz> <SNPs.vcf> <output.cov>
# specify output name without .gz; gzip is run at the end
# should work on either compressed or uncompressed inputs
# this could be repeated for multiple VCFs, e.g. from multiple SNP callers

# e.g.:
# sh filter_SNPs_from_covfile2.sh Dog4.CpG_report.merged_CpG_evidence.cov.gz Dog-4.g.vcf.gz Dog4.CpG_report.merged_CpG_evidence.SNPfilt.cov

# convert the (1-based) cov file to a 0-based BED format (subtract 1 from the second and third columns)
zcat $1 | awk -v OFS='\t' '$2=$2-1' | awk -v OFS='\t' '$3=$3-1' > tmp.bed
# run bedtools intersect with the -v option; this will remove any loci in the BED file that overlap with SNP positions in the VCF file
bedtools intersect -a tmp.bed -b $2 -v > tmp2.bed
# reformat the output BED file back into 1-based cov file (add 1 back second and third columns)
cat tmp2.bed | awk -v OFS='\t' '$2=$2+1' | awk -v OFS='\t' '$3=$3+1' > $3
# gzip the output (optional)
gzip $3

# with specific files:
# convert the (1-based) cov file to a 0-based BED format (subtract 1 from the second and third columns)
# zcat Dog4.CpG_report.merged_CpG_evidence.cov.gz | awk -v OFS='\t' '$2=$2-1' | awk -v OFS='\t' '$3=$3-1' > tmp.bed
# run bedtools intersect with the -v option; this will remove any loci in the BED file that overlap with SNP positions in the VCF file
# bedtools intersect -a tmp.bed -b Dog-4.g.vcf.gz -v > tmp2.bed
# reformat the output BED file back into 1-based cov file (add 1 back second and third columns)
# cat tmp2.bed | awk -v OFS='\t' '$2=$2+1' | awk -v OFS='\t' '$3=$3+1' > Dog4.CpG_report.merged_CpG_evidence.SNPfilt.cov
# gzip the output (optional)
# gzip Dog4.CpG_report.merged_CpG_evidence.SNPfilt.cov
