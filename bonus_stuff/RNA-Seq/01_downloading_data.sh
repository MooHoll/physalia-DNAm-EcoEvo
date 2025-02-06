#!/bin/bash

#SBATCH --job-name=downloading_data
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=14:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hollie_marshall@hotmail.co.uk
#SBATCH --account=evo-epi

### Include when testing
# SBATCH --partition=devel

### Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

# Load modules
module load sratoolkit/2.11.1

# Download data from NCBI
prefetch --option-file Accessions.txt

# Move all files out of the directories

# Make SRR files into fastq
for file in $(ls SRR*)
do
    fasterq-dump --split-files ${file}
done
