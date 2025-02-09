###---------------------------------------------------------------
# Examining genome-wide gene expression relationship with methylation
###---------------------------------------------------------------

library(readr)
library(reshape2)
library(Hmisc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcomp)

###---------------------------------------------------------------
# Read in the data

all_methylation_expression_data <- read_delim("all_methylation_expression_data.txt", 
                                       delim = "\t", escape_double = FALSE, trim_ws = TRUE)

###---------------------------------------------------------------
# Plot gene by gene scatter plot

ggplot(all_methylation_expression_data, aes(x=weighted_meth, y=logFPKM, colour = sex))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Take a look for just males
meth_exp_males <- all_methylation_expression_data[all_methylation_expression_data$sex=="male",]

ggplot(meth_exp_males, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#44AA99",size=2)+
  geom_smooth(method = "lm",size=2,colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Male")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)

# Take a look for just females
meth_exp_females <- all_methylation_expression_data[all_methylation_expression_data$sex=="female",]

ggplot(meth_exp_females, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#DDCC77",size=2)+
  geom_smooth(method = "lm",size=2, colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Female")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)


# Run some simple stats to see if methylation levels predict gene expression level
model1<-lm(FPKM ~ sex * weighted_meth, data=all_methylation_expression_data)
model2<-lm(FPKM ~ sex + weighted_meth, data=all_methylation_expression_data)
anova(model1,model2) # No interaction between sex and level of methylation
summary.lm(model2) # Sex not significant but meth does predict exp (just about)

###---------------------------------------------------------------
# The above could be the end of this analysis, the following graphs
# are representing the same data but in different ways and are a 
# personal preference
###---------------------------------------------------------------

###---------------------------------------------------------------
# Bin methylation into 100 bins and plot the mean of those bins

meth_exp_females$bins <- as.numeric(cut2(meth_exp_females$weighted_meth , g=100))
meth_exp_males$bins<-as.numeric(cut2(meth_exp_males$weighted_meth , g=100))

female<-as.data.frame(aggregate(meth_exp_females$FPKM, by=list(meth_exp_females$bins), mean))
female$logfpkm <- log10(female$x)
colnames(female)<-c("meth_bin","FPKM","logFPKM")
female$status <- "female"
plot(female$meth_bin~female$FPKM)

male<-as.data.frame(aggregate(meth_exp_males$FPKM, by=list(meth_exp_males$bins), mean))
male$logfpkm <- log10(male$x)
colnames(male)<-c("meth_bin","FPKM","logFPKM")
male$status <- "male"
plot(male$meth_bin~male$FPKM)

final_data <- merge(male, female, by="meth_bin")
final_data <- rbind(male, female)

ggplot(data=final_data, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point(size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))


###---------------------------------------------------------------
# Put methylation into even larger categories for a violin plot

# Define a function for confidence intervals
condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Make categories by methylation level and sex
all_methylation_expression_data$plot_cat <- paste(all_methylation_expression_data$meth_level,
                                                  all_methylation_expression_data$sex, sep="_")
ggplot(all_methylation_expression_data, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Methylation Category")+
  ylab("log(FPKM)")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("None_female","None_male","Low_female","Low_male","Medium_female",
                            "Medium_male","High_female","High_male"),
                   labels=c("None", "","Low","","Medium","","High",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


###---------------------------------------------------------------
# Make the same graphs but for autosomes and the X-chromosome

meth_exp_X <- all_methylation_expression_data[all_methylation_expression_data$chr=="DC3.0sc08",] #1428
meth_exp_A <- all_methylation_expression_data[!all_methylation_expression_data$chr=="DC3.0sc08",] #22546

# Scatters
meth_exp_males_X <- meth_exp_X[meth_exp_X$sex=="male",]
meth_exp_males_A <- meth_exp_A[meth_exp_A$sex=="male",]

ggplot(meth_exp_males_A, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#44AA99",size=2)+
  geom_smooth(method = "lm",size=2,colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Male: Autosomes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)

meth_exp_females_X <- meth_exp_X[meth_exp_X$sex=="female",]
meth_exp_females_A <- meth_exp_A[meth_exp_A$sex=="female",]

ggplot(meth_exp_females_A, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#DDCC77",size=2)+
  geom_smooth(method = "lm",size=2, colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Female: Autosomes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)


# Stats
model1<-lm(FPKM~sex*weighted_meth, data=meth_exp_X)
model2<-lm(FPKM~sex+weighted_meth, data=meth_exp_X)
anova(model1,model2) # No interaction
summary.lm(model2) # Nothing signif

# Stats
model1<-lm(FPKM~sex*weighted_meth, data=meth_exp_A)
model2<-lm(FPKM~sex+weighted_meth, data=meth_exp_A)
anova(model1,model2) # No interaction
summary.lm(model2) # Sex not sig but meth exp is


# Binned line graph
meth_exp_females_X$bins <- as.numeric(cut2(meth_exp_females_X$weighted_meth , g=100))
meth_exp_males_X$bins<-as.numeric(cut2(meth_exp_males_X$weighted_meth , g=100))

meth_exp_females_A$bins <- as.numeric(cut2(meth_exp_females_A$weighted_meth , g=100))
meth_exp_males_A$bins<-as.numeric(cut2(meth_exp_males_A$weighted_meth , g=100))

female_X<-as.data.frame(aggregate(meth_exp_females_X$FPKM, by=list(meth_exp_females_X$bins), mean))
female_X$logfpkm <- log10(female_X$x)
colnames(female_X)<-c("meth_bin","FPKM","logFPKM")
female_X$status <- "female"

female_A<-as.data.frame(aggregate(meth_exp_females_A$FPKM, by=list(meth_exp_females_A$bins), mean))
female_A$logfpkm <- log10(female_A$x)
colnames(female_A)<-c("meth_bin","FPKM","logFPKM")
female_A$status <- "female"

male_X<-as.data.frame(aggregate(meth_exp_males_X$FPKM, by=list(meth_exp_males_X$bins), mean))
male_X$logfpkm <- log10(male_X$x)
colnames(male_X)<-c("meth_bin","FPKM","logFPKM")
male_X$status <- "male"

male_A<-as.data.frame(aggregate(meth_exp_males_A$FPKM, by=list(meth_exp_males_A$bins), mean))
male_A$logfpkm <- log10(male_A$x)
colnames(male_A)<-c("meth_bin","FPKM","logFPKM")
male_A$status <- "male"

final_data_X <- merge(male_X, female_X, by="meth_bin")
final_data_X <- rbind(male_X, female_X)

final_data_A <- merge(male_A, female_A, by="meth_bin")
final_data_A <- rbind(male_A, female_A)

ggplot(data=final_data_X, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point(size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  ggtitle("X-Chromosome")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                      values=c("#DDCC77","#44AA99"))

ggplot(data=final_data_A, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point(size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  ggtitle("Autosomes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                      values=c("#DDCC77","#44AA99"))

# Binned violin plots
meth_exp_X$plot_cat <- paste(meth_exp_X$meth_level, meth_exp_X$sex, sep="_")
meth_exp_A$plot_cat <- paste(meth_exp_A$meth_level, meth_exp_A$sex, sep="_")

ggplot(meth_exp_X, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Methylation Category")+
  ylab("log(FPKM)")+
  ggtitle("X-Chromosome")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("None_female","None_male","Low_female","Low_male","Medium_female",
                            "Medium_male","High_female","High_male"),
                   labels=c("None", "","Low","","Medium","","High",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


ggplot(meth_exp_A, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Methylation Category")+
  ylab("log(FPKM)")+
  ggtitle("Autosomes")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("None_female","None_male","Low_female","Low_male","Medium_female",
                            "Medium_male","High_female","High_male"),
                   labels=c("None", "","Low","","Medium","","High",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())
