#------------------------------------------
# Make files for a pca from the vcf
#------------------------------------------

# Run an interactive job
salloc --partition=devel --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --time=1:0:0 --mem=16g
srun --pty bash -i

# Load modules
module load plink2
module load vcftools

# Make plink files
plink2 --vcf all_samples_filtered.recode.vcf  --allow-extra-chr --out daphnia_ER
plink2 --vcf all_samples_filtered.recode.vcf --make-pgen --allow-extra-chr --out daphnia_ER
plink2 --vcf all_samples_filtered.recode.vcf --freq --allow-extra-chr --out daphnia_ER

# Make eigenvectors for pca
plink2 --vcf all_samples_filtered.recode.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--read-freq daphnia_ER.afreq \
--max-alleles 2 \
--make-bed --pca 6 --out daphnia_ER