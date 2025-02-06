#!/bin/bash

# getting estimate of conversion efficiency entails aligning only to the chloroplast genome

N_CORES=2
READS=Col0_subsample.trimmed.fastq.gz
PREFIX=Col0_align_chloroplast

# index the TAIR10 chloroplast genome
mkdir bisulfite_genome_chloroplast
cp Arabidopsis_thaliana.TAIR10.chloroplast.fa bisulfite_genome_chloroplast/
bismark_genome_preparation bisulfite_genome_chloroplast

# run bismark align
bismark --prefix $PREFIX -p $N_CORES --bam bisulfite_genome_chloroplast $READS