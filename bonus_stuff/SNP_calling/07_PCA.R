# ------------------------------------------------------------
# Making PCA for SNPs
# ------------------------------------------------------------

library(tidyverse)
library(readr)
library(ggrepel)
library(ggplot2)

setwd("~/STUDENTS/PROJECT_MATERIALS/SNP_calling/Step_7")

# ------------------------------------------------------------

# Read in data from PLINK2
pca <- read_table2("./Cheat_files/daphnia_ER.eigenvec", col_names = T)
eigenval <- scan("./Cheat_files/daphnia_ER.eigenval")

# Fix up dataframe
names(pca)[1] <- "Line"
pca <- pca[,-2]
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Make percentage of variance explained by each PC
pve <- data.frame(PC = 1:6, pve = eigenval/sum(eigenval)*100)

# ------------------------------------------------------------

# Scree plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("Percentage variance explained") + theme_light()

# Check what the percentage variance explained adds to
cumsum(pve$pve)

# Make a vector of percentages for the axis label
percentage <- round(pve[,2],1)

# Accidental extra sample - remove
pca <- pca [-2,]

pca$Population <- c(rep("Experienced",3), rep("Naive",3))
pca$genotype <- c("6_2","7_5","9_20","36_01","53_01","54_01")

# ------------------------------------------------------------

# Plot PCA, lab line is super different as expected
ggplot(pca, aes(PC1, PC2, label = genotype, colour=Population)) + 
  geom_point(size = 6)+
  coord_equal() +
  geom_text(hjust="inward", vjust=2) +
  xlab(paste0("PC1 ",percentage[1],"%")) +
  ylab(paste0("PC2 ",percentage[2],"%")) +
  theme_bw()
