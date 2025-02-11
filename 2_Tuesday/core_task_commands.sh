# ---------------------------------------------------------
# Commands for the Core Tasks of short-read processing
# ---------------------------------------------------------

# Change directory to where you .pem file is located
# Log in using the ssh command
ssh -i <your_number>.pem user<number>@<IP_code>

# Activate the conda environment
conda activate /home/ubuntu/miniconda3/envs/short_read

# ---------------------------------------------------------
# 1a. Quality checking and adapter trimming
# ---------------------------------------------------------

# Get data into your area
cp ~/Share/1_Monday/raw_data_for_alignment/Col-0/*gz ./

# Subsample your data to a smaller file, this just makes it easier for us to work with
# This uses the BBtools package (already installed in the conda environment)
reformat.sh in=./SRR771698.fastq.gz out=Col0_subsample.fastq.gz samplerate=0.1

# Run the quality checking
# This uses the fastp package (already installed in the conda environment)
fastp --thread 2 --in1 Col0_subsample.fastq.gz --out1 Col0_subsample.trimmed.fastq.gz \
    -l 35 -h Col0_subsample_SE.fastp.html \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20

# Use filezilla or your terminal to get the HTML report file to look at on your computer
exit # to exit the server
cd ~/Desktop # somewhere convenient on your computer
scp -i <YOUR>.pem <user>@34.209.238.173://home/ubuntu/*.html ./

# Log back into the server
# Re-ativate the conda environment

# ---------------------------------------------------------
# 1b. Alignment to the reference genome
# ---------------------------------------------------------

# Get the genome folder (with the .fa genome file inside)
cp -r ~/Share/1_Monday/TAIR10_genome ./

# Prepare the genome for alignment
bismark_genome_preparation ./TAIR10_genome

# Align your file to the referece genome
bismark --multicore 2 --genome ./TAIR10_genome/ Col0_subsample.trimmed.fastq.gz

# Paired-end alignment for those working with the paired data
#bismark --multicore 2 --genome ./TAIR10_genome/ -1 SRR3301595_1.fastq.gz -2 SRR3301595_2.fastq.gz

# Examine the output report file
nano Col0_subsample.trimmed.fastq_bismark_bt2_PE_report.txt

# What is the alignment rate?
# What is the level of DNA methylation?

# ---------------------------------------------------------
# 2a. Conversion efficency:
# ---------------------------------------------------------

# Get and pepare the chloroplast genome
cp -r ~/Shared/2_Tuesday/TAIR10_chloroplast/ ./
bismark_genome_preparation ./TAIR10_chloroplast

# As above - single end alignment
bismark --prefix chloroplast --multicore 2 --genome ./TAIR10_genome/ Col0_subsample.trimmed.fastq.gz

# Paired end alignment
# bismark --prefix chloroplast --multicore 2 --genome ./TAIR10_chloroplast/ -1 SRR3301595_1.fastq.gz -2 SRR3301595_2.fastq.gz

# Examine the output files
nano chloroplast.Col0_subsample.trimmed.fastq_bismark_bt2_PE_report.txt 

# What is the alignment rate?
# What is the level of DNA methylation?

# ---------------------------------------------------------
# 2b. Methylation extraction
# ---------------------------------------------------------

# Deduplicate your .bam alignments (the ones aligned to the arabidopsis genome, not the chloroplast genome)
deduplicate_bismark Col0_subsample.trimmed.fastq_bismark_bt2.bam

# Extract the methylation data
bismark_methylation_extractor --multicore 3 --comprehensive --bedGraph --report --cytosine_report \
--genome_folder ./TAIR10_genome/ Col0_subsample.trimmed.fastq_bismark_bt2.bam

# Examine the M-bias plots on your computer using Filezilla or the below copy command
exit # to exit the server
cd ~/Desktop # somewhere convenient on your computer
scp -i <YOUR>.pem <user>@34.209.238.173://home/ubuntu/*.png ./ # edit this for your specific login details and path

# Log back into the server
# Re-ativate the conda environment

# Optional: merge strands
coverage2cytosine -o Col0 --merge_CpGs --genome_folder ./TAIR10_genome/ Col0_subsample.trimmed.fastq_bismark_bt2.bam

# ---------------------------------------------------------
# 2c. Optional: Accounting for SNPs
# ---------------------------------------------------------

# Copy the SNP files for your chosen strain from the /Shared/ area to your area
cp ~/Shared/2_Tuesday/SNP_calls/Col-0.g.vcf.gz ./

# Get the C/T and G/A SNPs we want to filter from a vcf
zcat Col-0.g.vcf.gz | grep -v "#" | awk '$4=="C"&&$5=="T,<NON_REF>"' | cut -f1,2 > Col-0_SNPs_CtoT.txt
zcat Col-0.g.vcf.gz | grep -v "#" | awk '$4=="G"&&$5=="A,<NON_REF>"' | cut -f1,2 > Col-0_SNPs_GtoA.txt

# Run the filtering script provided - coverage files end in .cov or .cov.gz
bash filter_SNPs_from_covfile.sh <your_coverage_file> Col-0_SNPs_CtoT.txt
bash filter_SNPs_from_covfile.sh <output_from_command_above> Col-0_SNPs_GtoA.txt
