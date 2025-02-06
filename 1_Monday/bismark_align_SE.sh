#!/bin/bash

# align SE read file to TAIR10 genome

N_CORES=2
INDEX=bisulfite_genome_TAIR10
READS=Col0_subsample.trimmed.fastq.gz
PREFIX=Col0_align_TAIR10

bismark --prefix $PREFIX -p $N_CORES --bam $INDEX $READS

# last two arguments are the folder containing the genome index and the read file to be aligned