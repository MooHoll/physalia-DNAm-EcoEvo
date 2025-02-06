### Snakemake workflows on ALICE2: Univeristy of Leicester HPC

Use a conda environment to call snakemake. 
Create a config file, Snakefile and .slm script.

conda create --name snakemake_env

conda activate snakemake_env

conda install -c bioconda snakemake

NOTE: If running tests need to run: snakemake --unlock when done, before submitting the job to the cluster,
otherwise it will kill it or put --nolock after snakemake when running

NOTE: also need to use bioconda to install all of the other software needed, like bismark