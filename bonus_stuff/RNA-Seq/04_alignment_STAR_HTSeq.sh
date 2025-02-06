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
#conda install bioconda::htseq

source ~/miniconda3/bin/activate alignment

### Genome prep
STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles ./phaw_5.0_Oxford.fa \
--sjdbGTFfile ./phaw_5.0_Oxford.gff3 \
--sjdbOverhang 99 \
--genomeChrBinNbits 15 \
--limitGenomeGenerateRAM=32000000000

### Alignment
for file in $(ls *_trim_1.fq.gz)
do
    base=$(basename $file "_trim_1.fq.gz")
    STAR \
    --runThreadN 24 \
    --genomeDir ./ \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}_trim_1.fq.gz ${base}_trim_2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${base}
done


### Counting reads
for file in $(ls *Aligned.sortedByCoord.out.bam)
do
    base=$(basename $file "Aligned.sortedByCoord.out.bam")
    htseq-count --format=bam  \
    ${base}Aligned.sortedByCoord.out.bam \
    phaw_5.0_Oxford.gff3 \
    > ${base}.counts
done