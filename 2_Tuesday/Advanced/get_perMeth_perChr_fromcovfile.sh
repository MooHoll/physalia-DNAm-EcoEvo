#!/bin/bash

# script to get average methylation % from each chromosome, taking input bismark.cov.gz file
# useful for checking BS conversion efficiency, e.g. via meth % of chloroplast
# run as:
# sh get_perMeth_perChr_fromcovfile.sh file.cov.gz

# extract chromosome names
zcat $1 | cut -f1 | sort -u > tmp1

# for each chromosome name
for i in `cat tmp1`;do
# extract rows from the cov file where the value of column 1 matches the current chromosome name,
# and calculate the average of values in column 4 (% methylation)
zcat $1 | awk -v i="$i" '$1==i' | awk '{ total += $4 } END { print total/NR }' >> tmp2
done
# print the chromosome names and average methylation %
paste tmp1 tmp2

# remove temporary files
rm tmp*