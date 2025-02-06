# In order to use Snakemake you need a conda environment
# Here is how to make one

# Create an environment
conda create --name snakemake_env

# Activate it
conda activate snakemake_env

# Install all the software you need
conda install -c bioconda snakemake
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda picard
conda install -c bioconda freebayes
conda install -c bioconda vcftools