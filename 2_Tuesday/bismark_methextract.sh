#!/bin/bash

BAM=Col0_subsample.deduplicated.bam

bismark_methylation_extractor --bedGraph --merge_non_CpG --comprehensive --gzip $BAM
# add --paired-end option if dealing with paired end data

# default bismark coverage file is for CpG context only
# to get coverage file only for CHH context, for example:
# bismark2bedGraph -o CHH_context_Col0_subsample CHH_context_Col0_align_TAIR10.Col0_subsample.trimmed_bismark_bt2.txt.gz