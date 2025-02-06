#-------------------------------------------------------------------
# Blast to identify DNMTs/TET in your species
#-------------------------------------------------------------------

# Get protein sequences of all genes of your species
conda create -n agat_env
conda activate agat_env
conda install -c bioconda agat

agat_sp_extract_sequences.pl --gff phaw_5.0_Oxford.gff3 -f phaw_5.0_Oxford.fa -p -o parhyale_proteins.fa

# Make a blast databases
module load blast-plus/2.13.0-5o3kbvq

makeblastdb -in parhyale_proteins.fa -parse_seqids -dbtype prot


# Using http://v2.insect-genome.com/Pcg, could download the protein .fa seq in one file
# for DNMT1 and DNMT3a for 321 and 110 species respectively 
# for TET >1000 sequences from >700 species
# (choose species specific reference DNMTS/TET, the above are for insects)
makeblastdb -in dnmt1.fa -parse_seqids -dbtype prot
makeblastdb -in dnmt3a.fa -parse_seqids -dbtype prot
makeblastdb -in tet.fa -parse_seqids -dbtype prot

# reciprocal for DNMT1
blastp -query parhyale_proteins.fa \
-db "dnmt1.fa" -evalue 1e-3 -max_target_seqs 5 \
-outfmt 6 -out parhyale_DNMT1.txt

blastp -query dnmt1.fa \
-db "parhyale_proteins.fa" -evalue 1e-3 -max_target_seqs 5 \
-outfmt 6 -out parhyale_DNMT1_reciprocol.txt

module load R/4.3.1
R
library(readr)
data <- read_delim("parhyale_DNMT1_reciprocol.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=F)
table(data$X2)

# reciprocal for DNMT3
blastp -query parhyale_proteins.fa \
-db "dnmt3a.fa" -evalue 1e-3 -max_target_seqs 5 \
-outfmt 6 -out parhyale_DNMT3.txt

blastp -query dnmt3a.fa \
-db "parhyale_proteins.fa" -evalue 1e-3 -max_target_seqs 5 \
-outfmt 6 -out parhyale_DNMT3_reciprocol.txt

R
library(readr)
data <- read_delim("parhyale_DNMT3_reciprocol.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=F)
table(data$X2)

# reciprocal for TET
blastp -query parhyale_proteins.fa \
-db "tet.fa" -evalue 1e-3 -max_target_seqs 5 \
-outfmt 6 -out parhyale_TET.txt

blastp -query tet.fa \
-db "parhyale_proteins.fa" -evalue 1e-3 -max_target_seqs 5 \
-outfmt 6 -out parhyale_TET_reciprocol.txt

R
library(readr)
data <- read_delim("parhyale_TET_reciprocol.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=F)
table(data$X2)