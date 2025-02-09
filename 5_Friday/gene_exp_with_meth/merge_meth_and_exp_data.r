###-----------------------------------------------------------------
# Make one dataframe with gene expression and methylation together
###-----------------------------------------------------------------

library(readr)
library(reshape2)
library(dplyr)
library(tidyr)

###-----------------------------------------------------------------
# Merge the weighted methylation data with the differential gene methylation

weighted_meth <- read_delim("Dcitri_weighted_meth_genes_only.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

diff_meth_genes <- read_delim("Dcitri_differentially_methylated_genes.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
diff_meth_genes <- diff_meth_genes[,c(1,3,10)]

methylation <- merge(weighted_meth, diff_meth_genes, by=c("chr","gene_id"), all=T)
methylation$hypermethylated[is.na(methylation$hypermethylated)] <- "none"

methylation$male_all <- paste(methylation$male, methylation$male_category, sep="_")
methylation$female_all <- paste(methylation$female, methylation$female_category, sep="_")

methylation <- methylation[,-c(3:6)]
methylation <- reshape2::melt(methylation, id.vars=c("chr","gene_id","hypermethylated"))
methylation <- methylation[!duplicated(methylation),]
methylation <- separate(methylation, value, into = c("weighted_meth","meth_level"), sep="_")
colnames(methylation)[4] <- "sex"
methylation$sex <- gsub("_all","", methylation$sex)

###-----------------------------------------------------------------
# Read in gene expression data

exp_data <- read_delim("Dcitri_gene_expression.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
exp_data <- exp_data[,-c(2:7,11:15)]
exp_data <- reshape2::melt(exp_data, id.vars=c("gene_id","diff_exp","category", "log2FoldChange"))
colnames(exp_data) <- c("gene_id","diff_exp","category","logFC","sex","FPKM")
exp_data$sex <- gsub("_fpkm_mean","", exp_data$sex)

###-----------------------------------------------------------------
# Bring both data sets together

meth_exp <- merge(exp_data, methylation, by=c("gene_id","sex"), all=T) #35924 genes
meth_exp <- meth_exp[!is.na(meth_exp$weighted_meth) & !is.na(meth_exp$FPKM),] #23964 genes

# Avoid complete separation in later models
meth_exp$FPKM[meth_exp$FPKM==0] <- 1
meth_exp$logFPKM <- log(meth_exp$FPKM)
meth_exp$logFPKM[meth_exp$logFPKM == -Inf] <- 0

meth_exp$weighted_meth <- as.numeric(meth_exp$weighted_meth)

write.table(meth_exp, file="all_methylation_expression_data.txt",
          col.names = T, row.names = F, quote = F, sep = "\t")