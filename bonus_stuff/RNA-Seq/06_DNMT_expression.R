#---------------------------------------------------------
# Expression levels of DNMT genes
#---------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Differential_expression")

library(readr)
library(reshape2)
library(ggplot2)
#---------------------------------------------------------

# This file comes from the DESeq2 script
all_gene_expression_data <- read_delim("all_gene_expression_data.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)

# DNMT1 -1 Dcitr08g10610.1
# DNMT1 -2 Dcitr08g05090.1
# DNMT2 Dcitr07g04270.1

head(all_gene_expression_data)

dnmt1_1 <- all_gene_expression_data[all_gene_expression_data$gene_id=="Dcitr08g10610.1",]
# Unbiased

dnmt1_2 <- all_gene_expression_data[all_gene_expression_data$gene_id=="Dcitr08g05090.1",]
# Not present

dnmt2 <- all_gene_expression_data[all_gene_expression_data$gene_id=="Dcitr07g04270.1",]
# Unbiased

dnmt1_1 <- dnmt1_1[,c(1:7)]
dnmt1_1 <- melt(dnmt1_1)
dnmt1_1$Sex <- c("Female","Female","Female","Male","Male","Male")

ggplot(dnmt1_1, aes(x=Sex, y=value, fill=Sex))+
  geom_boxplot()+
  geom_jitter(size=5)+
  theme_bw()+
  xlab("Sex")+
  ylab("Expression Level (FPKM)")+
  ggtitle("DNMT1")+
  scale_fill_manual(limits=c("Female", "Male"),
                    values = c("#DDCC77","#44AA99"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))


dnmt2 <- dnmt2[,c(1:7)]
dnmt2 <- melt(dnmt2)
dnmt2$Sex <- c("Female","Female","Female","Male","Male","Male")

ggplot(dnmt2, aes(x=Sex, y=value, fill=Sex))+
  geom_boxplot()+
  geom_jitter(size=5)+
  theme_bw()+
  xlab("Sex")+
  ylab("Expression Level (FPKM)")+
  ggtitle("DNMT2")+
  scale_fill_manual(limits=c("Female", "Male"),
                    values = c("#DDCC77","#44AA99"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))

# Also take a look at the dsx and fru genes
head(all_gene_expression_data)

dsx <- all_gene_expression_data[all_gene_expression_data$gene_id=="Dcitr03g16970.1",]
# Not present

fru <- all_gene_expression_data[all_gene_expression_data$gene_id=="Dcitr01g04580.1",]
# Unbiased

fru <- fru[,c(1:7)]
fru <- melt(fru)
fru$Sex <- c("Female","Female","Female","Male","Male","Male")

ggplot(fru, aes(x=Sex, y=value, fill=Sex))+
  geom_boxplot()+
  geom_jitter(size=5)+
  theme_bw()+
  xlab("Sex")+
  ylab("Expression Level (FPKM)")+
  ggtitle("Fruitless")+
  scale_fill_manual(limits=c("Female", "Male"),
                    values = c("#DDCC77","#44AA99"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))
