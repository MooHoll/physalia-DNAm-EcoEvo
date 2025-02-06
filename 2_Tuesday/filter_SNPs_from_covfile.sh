#!/bin/bash

# script to generate SNP-filtered bismark coverage file, given (1) the cov.gz file and (2) a list of SNP positions in tab-separated format (chr, loc)
# run as:
# filter_SNPs_from_covfile.sh bismark.cov.gz SNPs.txt

awk '{print $1"_"$2}' $2 > SNPIDs.tmp
# make new .cov file with added site IDs (in format 'chr_loc')
zcat $1 | awk '{print $1"_"$2, $0}' > tmp1
# remove lines from the new file with site ID matching SNP ID
awk 'NR == FNR {a[$1]; next} !($1 in a)' SNPIDs.tmp tmp1 > tmp2
# remove site ID column to convert back to original file format; command to the right of the pipe converts the file back to tab-separated format
awk '{print $2,$3,$4,$5,$6,$7}' tmp2 | awk -v OFS='\t' '{ $1=$1; print }' > ${1}.SNPfilt.cov
gzip ${1}.SNPfilt.cov
# remove site ID file
rm *tmp*