#!/bin/bash

# index the TAIR10 genome
mkdir bisulfite_genome_TAIR10
cp Arabidopsis_thaliana.TAIR10.dna.toplevel.fa bisulfite_genome_TAIR10/
bismark_genome_preparation bisulfite_genome_TAIR10