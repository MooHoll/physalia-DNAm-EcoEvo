#!/bin/bash

# remove PCR duplicates from BAM file

BAM=Col0_align_TAIR10.Col0_subsample.trimmed_bismark_bt2.bam
PREFIX=Col0_subsample

deduplicate_bismark --bam --outfile $PREFIX $BAM
# add --paired if data were paired-end

# Note that this may not be perfect, given that quality trimming may affect some duplicates differently