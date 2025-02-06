#!/bin/bash

#SBATCH --job-name=snakemake
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

source ~/miniconda3/bin/activate snakemake_env

# keep going means if there is an error it will run independend jobs still
# no lock needed to don't need to keep unlocking the working directory
snakemake --rerun-incomplete --keep-going --nolock