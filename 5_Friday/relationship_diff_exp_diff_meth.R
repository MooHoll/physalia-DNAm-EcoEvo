# ---------------------------------------------------------------
# Relationship between differential methylation and differential exp
# ---------------------------------------------------------------

library(readr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(tidyr)
library(multcomp)
library(UpSetR)
library(stringr)

# ---------------------------------------------------------------
# Read in all data

all_methylation_expression_data <- read_delim("all_methylation_expression_data.txt", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

# Calculate the differenece in weighted methylation between the sexes
weighted_meth <- all_methylation_expression_data[,c(1,2,9)]
weighted_meth <- weighted_meth[!duplicated(weighted_meth),]
weighted_meth <- spread(weighted_meth, sex, weighted_meth)
weighted_meth$meth_diff <- as.numeric(weighted_meth$male - weighted_meth$female)
weighted_meth <- weighted_meth[,c(1,4)]

all <- merge(all_methylation_expression_data, weighted_meth, by = "gene_id", all=T)

# ---------------------------------------------------------------
# Are there any genes which are both differentially expressed and differentially methylated?

genes <- all[!(all$diff_exp=="no") & !(all$hypermeth_cat=="none"),] # 0


# ---------------------------------------------------------------
# Scatter plots of all differential data coloured by differential gene expression

# Remove one non-significant outlier to make a more visually appealing graph
all$logFC <- as.numeric(all$logFC)
all <- all[!(all$logFC < -15),]

# Scatter plot all data and colour differentially expressed genes
ggplot(all, aes(x=meth_diff, y=logFC, colour=diff_exp))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  theme_bw()+
  ggtitle("Differentially Expressed Genes")+
  scale_colour_manual("", breaks=c("no", "yes"),
                      values = c("black","red"),
                      labels= c("No", "Yes"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)

# ---------------------------------------------------------------
# Scatter plots of all differential data coloured by differential methylation

# Add categories for easier plot colouring
all$diff_meth_cat <- "yes"
all$diff_meth_cat[all$hypermethylated=="none"] <- "no"

# Order so we can see the coloured dots on the graph
ordered <- all[order(all$diff_meth_cat),]

ggplot(ordered, aes(x=meth_diff, y=logFC, colour=diff_meth_cat))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  theme_bw()+
  ggtitle("Differentially Methylated Genes")+
  scale_colour_manual("", breaks=c("no", "yes"),
                      values = c("black","red"),
                      labels= c("No", "Yes"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)



# ---------------------------------------------------------------
# Look at the same data in violin plot format
# ---------------------------------------------------------------

# Define confidence interval function
condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# ---------------------------------------------------------------
# Examine the expression levels of differentially methylated genes

all$plot_cat <- paste(all$hypermethylated, all$sex, sep="_")

ggplot(all, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Differential Methylation Category")+
  ylab("log(FPKM)")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("none_female","none_male","female_female","female_male",
                            "male_female","male_male"),
                   labels=c("None", "None","Female Hypermethylated \nExon","", "Male Hypermethylated \nExon",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=16, hjust=0.01),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


# Are there significant differences in gene expression between genes with
# differentially methylated exons between the sexes?
all$sex <- as.factor(all$sex)
all$hypermethylated <- as.factor(all$hypermethylated)

model1<-lm(FPKM ~ sex * hypermethylated, data=all)
model2<-lm(FPKM ~ sex + hypermethylated, data=all)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Nope

# ---------------------------------------------------------------
# Examine the methylation levels of differentially expressed genes

all$log_meth <- log(all$weighted_meth)

# Make larger categories of types of differentially expressed genes for nicer plots
all$higher_cat <- "unbiased"
all$higher_cat[all$category=="male_biased"] <- "male_biased"
all$higher_cat[all$category=="male_limited"] <- "male_biased"
all$higher_cat[all$category=="male_biased_extreme"] <- "male_biased"
all$higher_cat[all$category=="female_biased"] <- "female_biased"
all$higher_cat[all$category=="female_limited"] <- "female_biased"
all$plot_column <- paste0(all$sex,"_",all$higher_cat)

ggplot(all, aes(x=plot_column, y=log_meth)) + 
  geom_violin( aes(fill=sex))+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  xlab("Differential Expression Category")+
  ylab("log(Weighted Methylation)")+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("female_unbiased","male_unbiased", "female_female_biased",
                            "male_female_biased","female_male_biased","male_male_biased"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Stats
model1<-lm(weighted_meth ~ sex * higher_cat, data=all)
model2<-lm(weighted_meth ~ sex + higher_cat, data=all)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Sig: differentially expressed genes have lower methylation than non-differentially expressed genes
