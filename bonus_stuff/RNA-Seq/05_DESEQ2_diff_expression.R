#------------------------------------------------------------------
# Differential Gene Expression
#------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Differential_expression")
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(apeglm)
library(tximportData)
library(tximport)

#------------------------------------------------------------------
# Make sample metadata
sample <- c("F1","F2","F3", "M1","M2","M3")
sex <- c("Female","Female","Female","Male","Male","Male")
file_name <- list.files("../Counts/", pattern="*.genes.results")
sample_info <- data.frame(sample, sex, file_name)
rownames(sample_info) <- sample_info$sample

# Read in files using tximport
files <- file.path("../Counts", paste0(sample_info$file_name))
names(files) <- paste0(sample_info$sample)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
txi.rsem$length[txi.rsem$length == 0] <- 1

# Make deseq2 object
dds <- DESeqDataSetFromTximport(txi.rsem,
                                   colData = sample_info,
                                   design = ~ sex)

#------------------------------------------------------------------
# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) # 12420/19049

# rlog transform counts (by average Tx length and correcting for library size)
rld = rlog(dds, blind=FALSE)

#------------------------------------------------------------------
# PCA plot
data = plotPCA(rld, intgroup = c("sex"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=sex)) + geom_point(size=14) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed(ratio=4,clip = "on")+
  geom_text_repel(aes(label=name), size=10,show.legend=FALSE, 
                  point.padding = 2, box.padding = 1,
                  segment.color = 'transparent') +
  scale_colour_manual("", breaks=c("Female","Male"),
                      values = c("#DDCC77","#44AA99"))+
  theme_bw()+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24),
        legend.text=element_text(size=24))

# Make a scree plot
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))] # Top variable genes
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- 100*(pca$sdev^2 / sum( pca$sdev^2 )) # the contribution to the total variance for each component

scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:6)
colnames(scree_plot)<-c("Variance","PC Number")

ggplot(scree_plot, mapping=aes(x=`PC Number`, y=Variance))+
  geom_bar(stat="identity")+
  theme_bw()+
  ylab("% Varience")+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24),
        legend.text=element_text(size=24))

# Heatmap of the top variable genes
my_colors = list(
  sex = c(Female = "#DDCC77", Male ="#44AA99"))
mat = assay(rld)[ select, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[,c("sex"),drop=FALSE])
colnames(df2)<-c("sex")

pheatmap(mat, annotation_col=df2,
         show_rownames = F,
         fontsize = 16,
         annotation_colors = my_colors)
# F1 weird genes: Dcitr08g09350.1, Dcitr06g09040.1, Dcitr06g09050.1, Dcitr01g10760.1, Dcitr08g09360.1, Dcitr08g09320.1, Dcitr08g09380.1


# two 1st samples plotted against each other to check consistency (for rlog and log2) 
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(1,2)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(1,2)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )

# estimate size factors = normalize for dispersion
dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

# check sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 18)

# check sample distances using the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         fontsize = 18)

#------------------------------------------------------------------
# differential expression (used Benjamini-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)
resultsNames(dds)

# Male = +Log2FC and Female = -Log2FC
res=results(dds, name="sex_Male_vs_Female")
summary(res)

#distribution of coefficents of the model
plotMA(res, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
       cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log2 Fold Change')

# plot of p-vals excluding genes with very small counts
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, main="")

# Order the results by fold change to see the biggest changing genes
res_ordered=res[order(res$log2FoldChange),]
head(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered$gene<-rownames(res_ordered)
rownames(res_ordered)<-c()
write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_all_genes.txt",
            sep="\t", quote = F, col.names = T, row.names = F)

#out of 12420 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 1.5 | log2FoldChange < -1.5)
res_significant <- subset(res_significant, padj < 0.05)
nrow(res_significant) #1259 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 1.5,]) #1164 upregulated in males
nrow(res_significant[res_significant$log2FoldChange < -1.5,]) #95 upregulated in females

#------------------------------------------------------------------
# heatmap of top differentially expressed
n=50
topdiff = head(c(1:nrow(res))[order(res$padj)],n)

my_colors = list(
  sex = c(Female = "#DDCC77", Male ="#44AA99"))


mat = assay(rld)[ topdiff, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[,c("sex"),drop=FALSE])
colnames(df2)<-c("sex")

pheatmap(mat, annotation_col=df2,
         show_rownames = F,
         fontsize = 16,
         annotation_colors = my_colors)

#------------------------------------------------------------------
# Plot diff expressed
n=12
selGenes = head(rownames(res)[order(res$padj)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("sex"), returnData=TRUE))))
ggplot(data, aes(x=sex, y=count, fill=sex)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Sex") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))+
  scale_fill_manual("", breaks=c("Female","Male"),
                    values = c("#DDCC77","#44AA99"))

#------------------------------------------------------------------
# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- row.names(res_df)

res_significant_df <- as.data.frame(res_significant)
res_significant_df$gene <- row.names(res_significant_df)

res_df$sig <- "no"
res_df$sig[res_df$gene %in% res_significant_df$gene] <- "yes"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-10,10)+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none")

#------------------------------------------------------------------
# Pull out a list of just the differentially expressed genes for the supplementary



