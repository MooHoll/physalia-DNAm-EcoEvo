#!/bin/bash

#PBS -N gatk_alternate_reference
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load gatk/3.8
module load java/1.8

# Remember the genome needs a .fai and .dict
# For .dict
#module load picard/2.6.0 
#java -Xms128m -Xmx1024m -Djava.io.tmpdir=$TMPDIR \
#-jar /cm/shared/apps/picard/2.6.0/picard.jar \
#CreateSequenceDictionary \
#R=PGA_assembly_DmagnaV3.fasta \
#O=PGA_assembly_DmagnaV3.dict

# Define file paths
REF_FILE=/scratch/monoallelic/hm257/hm343_stuff/genome/PGA_assembly_DmagnaV3.fasta                                                                             
GATK=/cm/shared/apps/gatk/3.8/GenomeAnalysisTK.jar

# create the directory where the output files are to be written   
OUTPUT=alternate_references                                                                                                                                   
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Run GATK to make the new alternative reference genomes
for file in $(ls *recode.vcf)
do
    base=$(basename ${file} ".recode.vcf")
    java -jar ${GATK}\
    -T FastaAlternateReferenceMaker \
    -R ${REF_FILE} \
    -V ${file} \
    -o ${OUTPUT}/${base}.fasta 
done