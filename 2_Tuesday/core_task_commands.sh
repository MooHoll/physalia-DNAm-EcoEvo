# ---------------------------------------------------------
# Commands for the Core Tasks: Day 2
# ---------------------------------------------------------

# ---------------------------------------------------------
# 1. Conversion efficency:
# ---------------------------------------------------------

# Activate the conda environment
conda activate /home/ubuntu/miniconda3/envs/short_read

# Prepare the chloroplast genome
cp -r ~/Shared/2_Tuesday/TAIR10_chloroplast/ ./
bismark_genome_preparation ./TAIR10_chloroplast

# Paired end alignment
bismark --multicore 3 --genome ./TAIR10_chloroplast/ -1 SRR3301595_1.fastq.gz -2 SRR3301595_2.fastq.gz
# Single end alignment
bismark --multicore 3 --genome ./TAIR10_genome/ SRR771746.fastq.gz 

# Examine the output files
nano chloroplast.SRR4295457_1.fastq_bismark_bt2_PE_report.txt 

# ---------------------------------------------------------
# 1. Methylation extraction
# ---------------------------------------------------------

# Deduplicate your .bam alignments (the ones aligned to the arabidopsis genome, not the chloroplast genome)
deduplicate_bismark SRR771746_bismark_bt2.bam

# Extract the methylation data
bismark_methylation_extractor --multicore 3 --comprehensive --bedGraph --report --cytosine_report \
--genome_folder TAIR10_genome/ SRR771746_bismark_bt2.bam

# Examine the M-bias plots
exit # to exit the server
cd ~/Desktop # somewhere convenient on your computer
scp -r -i <YOUR>.pem <user>@34.209.238.173://home/your/path/to/the/files/*.png ./ # edit this for your specific login details and path

# Log back into the server

# Optional: merge strands
coverage2cytosine -o SRR771746 --merge_CpGs --genome_folder TAIR10_genome/ SRR771746_bismark_bt2.bam

# ---------------------------------------------------------
# 3. Accounting for SNPs
# ---------------------------------------------------------

# Copy the SNP files for your chosen strain from the /Shared/ area to your area

# Get the C/T and G/A SNPs we want to filter from a vcf
zcat Col-0.g.vcf.gz | grep -v "#" | awk '$4=="C"&&$5=="T,<NON_REF>"' | cut -f1,2 > Col-0_SNPs_CtoT.txt
zcat Col-0.g.vcf.gz | grep -v "#" | awk '$4=="G"&&$5=="A,<NON_REF>"' | cut -f1,2 > Col-0_SNPs_GtoA.txt

# Run the filtering script provided
bash filter_SNPs_from_covfile.sh <your_coverage_file> Col-0_SNPs_CtoT.txt
bash filter_SNPs_from_covfile.sh <output_from_command_above> Col-0_SNPs_GtoA.txt