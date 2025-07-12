## -------------------------------------------------------------------------
## Pesticide vs Pristine
## -------------------------------------------------------------------------
setwd("~/PROJECTS/Daphnia_ETT_ER/1_Epigenetics_Through_Time/2023_analyses/WGBS/2_merged_cpgs")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Get a methylkit object for all samples
sample.list <- list("pesticide_6_2.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_6_3.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_7_3.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_7_5.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_7_5_4.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_8_5_3.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_9_20.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_9_5_1.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_9_5_3.CpG_report.merged_CpG_evidence.cov",
                    "pesticide_9_6.CpG_report.merged_CpG_evidence.cov",
                    "pristine_36_01.CpG_report.merged_CpG_evidence.cov",
                    "pristine_36_02.CpG_report.merged_CpG_evidence.cov",
                    "pristine_48_01.CpG_report.merged_CpG_evidence.cov",
                    "pristine_48_02.CpG_report.merged_CpG_evidence.cov",
                    "pristine_53_01.CpG_report.merged_CpG_evidence.cov",
                    "pristine_54_01.CpG_report.merged_CpG_evidence.cov",
                    "pristine_54_02.CpG_report.merged_CpG_evidence.cov",
                    "pristine_74_01.CpG_report.merged_CpG_evidence.cov",
                    "pristine_77_01.CpG_report.merged_CpG_evidence.cov",
                    "pristine_88_01.CpG_report.merged_CpG_evidence.cov")


CPGRaw <- methRead(sample.list, 
                   sample.id = list("pesticide_6_2",
                                    "pesticide_6_3",
                                    "pesticide_7_3",
                                    "pesticide_7_5",
                                    "pesticide_7_5_4",
                                    "pesticide_8_5_3",
                                    "pesticide_9_20",
                                    "pesticide_9_5_1",
                                    "pesticide_9_5_3",
                                    "pesticide_9_6",
                                    "pristine_36_01",
                                    "pristine_36_02",
                                    "pristine_48_01",
                                    "pristine_48_02",
                                    "pristine_53_01",
                                    "pristine_54_01",
                                    "pristine_54_02",
                                    "pristine_74_01",
                                    "pristine_77_01",
                                    "pristine_88_01"),
                   assembly="D.magna",
                   treatment=c(rep(0,9), rep(1,10)),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= path)

# Filter by coverage
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

# Select only CpGs found in all samples
meth_all_data <- unite(filtered_data, destrand=F) 
nrow(meth_all_data) 
# 1,021,209

# Pull data into a useable dataframe
df_meth_all <- getData(meth_all_data)

# remove first few uninformative columns
df_meth_all <- df_meth_all[,-c(1:4)]

# subset data into individual dataframes per sample
all_data <- lapply(seq(1, ncol(df_meth_all), by=3), function(i) 
  df_meth_all[i: pmin((i+1), ncol(df_meth_all))])

# Define binomial test
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}

# Set up output dataframe
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

# Run binomial 
for (i in seq_along(all_data)) {
  all_data[[i]] <- all_data[[i]][complete.cases(all_data[[i]]),]
  colnames(all_data[[i]]) <- c("CT", "Ccount")
  all_data[[i]]$pVal <- mapply(bt, all_data[[i]]$Ccount, all_data[[i]]$CT)
  all_data[[i]]$FDR <- p.adjust(all_data[[i]]$pVal, method = "BH", n = length(all_data[[i]]$pVal))
  all_data[[i]]$row <- as.numeric(rownames(all_data[[i]]))
  dfmeth <- subset(all_data[[i]], all_data[[i]]$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

# Get positions which are methylated in at least one sample
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 81282

# Keep only these positions for diff meth test
subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# Write out data for later use
subset_methBase_out <- getData(subset_methBase)
write.table(subset_methBase_out, file="sites_for_PCA.txt", sep="\t", quote=F, row.names = F, col.names = T)


## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(subset_methBase, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

PCA_data1$Condition <- c(rep("Experienced", 9), rep("Naive", 10))

percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )

#pdf(file = "./PCA_Bain_PoO_female.pdf",width=8.5, height =6)
ggplot(PCA_data1, aes(PC1, PC2, colour=Condition))+
  geom_point(size=14)+
  theme_bw()+
  xlab(paste0("PC1:",percentage[1],"variance")) +
  ylab(paste0("PC2:",percentage[2],"variance")) +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        legend.text=element_text(size=30),
        legend.title=element_blank(),
        title=element_text(size=26))+
 # ggtitle("PoO vs Bain")+
  scale_colour_manual(breaks=c("Experienced","Naïve"),
                      values=c("purple","yellow2"))
#dev.off()

# Scree plot
#PCASamples(subset_methBase, screeplot=T, obj.return = T)

## -------------------------------------------------------------------------
# Differential methylation of CpGs

diff_meth <- calculateDiffMeth(subset_methBase, method='qvalue')
diff_meth_filtered <- getMethylDiff(diff_meth, difference=15, qvalue=0.01)
nrow(diff_meth_filtered)
# 819

head(diff_meth_filtered)

write.csv(diff_meth_filtered, file="blah.csv")