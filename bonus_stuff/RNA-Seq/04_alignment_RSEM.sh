#!/bin/bash

#SBATCH --job-name=alignment
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=00:05:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hollie_marshall@hotmail.co.uk
#SBATCH --account=evo-epi

### Include when testing
# SBATCH --partition=devel

### Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

### Make a conda environment
#conda create -n alignment
#conda activate alignment
#conda install bioconda::star
#conda install bioconda::samtools
#conda install bioconda::rsem
#conda install bioconda::cufflinks

source ~/miniconda3/bin/activate alignment

# Convert gff to gtf for RSEM
gffread Daphnia_magna_LRV0_1.gff3 -T -o Daphnia_magna_LRV0_1.gtf

# Prepapre genome
rsem-prepare-reference --gtf Daphnia_magna_LRV0_1.gtf --star -p 1 Daphnia_magna_LRV0_1.scaffolds.fa daphnia

# Run alignment and read counting
for file in $(ls *_1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    rsem-calculate-expression --star --star-gzipped-read-file -p 16 --star-output-genome-bam \
    --paired-end ${base}_1.fq.gz ${base}_2.fq.gz \
    ./daphnia ${base}
done
