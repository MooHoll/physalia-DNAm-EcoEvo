#!/bin/bash

#SBATCH --job-name=whatever
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
module load py-multiqc/1.14-ateq76h
module load fastqc/0.12.1-hkgpcde

# Run fastqc to look at the quality of the data
for file in $(ls *.fq.gz)
do
	fastqc -t 24 ${file}
done

# Run multiqc to make a nice report
multiqc ./