#!/bin/bash

# align PE read files to TAIR10 genome

N_CORES=2
INDEX=bisulfite_genome_TAIR10
READ1=Dog4_subsample_1.trimmed.fastq.gz
READ2=Dog4_subsample_2.trimmed.fastq.gz
PREFIX=Dog4_align_TAIR10

bismark --prefix $PREFIX -p $N_CORES --bam $INDEX -1 $READ1 -2 $READ2

# last three arguments are the folder containing the genome index and the two paired read files to be aligned