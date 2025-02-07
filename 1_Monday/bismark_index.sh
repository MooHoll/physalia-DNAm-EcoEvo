#!/bin/bash

# index the TAIR10 genome
mkdir bisulfite_genome_TAIR10
cp ~/Share/1_Monday/TAIR10_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa bisulfite_genome_TAIR10
bismark_genome_preparation bisulfite_genome_TAIR10
