## -------------------------------------------------------------------------
# Weighted methylation per genomic feature of interest
## -------------------------------------------------------------------------

# This is quite an intensive script and may need running on a cluster for your real samples
library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)

## -------------------------------------------------------------------------

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*small_version.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("sample1","sample2","sample3")
names(samples) <- sample_names

# Read in features with start/end
annotation <- read_table2("random_species_annotations_small_version.bed")

# Tell R how many cores it has available to use
registerDoParallel(cores = 2)

# Calculate weighted meth for each feature for each sample (Windows)
foreach(i = seq_along(samples), .packages=c("sqldf","doBy","dplyr"),
        .export = ls(globalenv())) %dopar% {
  df <- samples[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id
                    FROM df AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$gene_id),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + feature + gene_id + start + end, data=output, FUN=sum) 
  check$weightedMeth <- (check$count_c.sum)/(check$total_coverage.sum)
  myfile <- file.path("./", paste0(names(samples[i]),"_","weighted_meth.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}

# Calculate weighted meth for each feature for each sample (Mac/Linux)
#foreach(i = seq_along(samples))) %dopar% {
#          df <- samples[[i]]
#          df <- subset(df, total_coverage > 5)
#          output <- sqldf("SELECT sample.chr,
#                    sample.cpg,
#                    sample.count_c,
#                    sample.total_coverage,
#                    annot.chr,
#                    annot.feature,
#                    annot.start,
#                    annot.end,
#                    annot.gene_id
#                    FROM df AS sample
#                    LEFT JOIN annotation AS annot
#                    ON sample.chr = annot.chr
#                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
#          output <- output[!is.na(output$gene_id),]
#          output <- output[,-c(1,2)]
#          check <- summaryBy(total_coverage + count_c ~ chr + feature + gene_id + start + end, data=output, FUN=sum) 
#          check$weightedMeth <- (check$count_c.sum)/(check$total_coverage.sum)
#          myfile <- file.path("./", paste0(names(samples[i]),"_","weighted_meth.txt"))
#          write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
#        }

## -------------------------------------------------------------------------

# Make one dataframe with all of the samples together
file.list = list.files(("./"),pattern="*weighted_meth.txt")
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}
samples_wm <- lapply(file.list, read_file1)
sample_names <- list("sample1","sample2","sample3")
names(samples_wm) <- sample_names

# Add sample name as a column
samples_with_name <- lapply(names(samples_wm),
                        function(current_name)transform(samples_wm[[current_name]],
                        sample = current_name))

# Merge dataframes
all <- as.data.frame(bind_rows(samples_with_name))

write.table(all, file="weighted_meth_all_samples.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")


## -------------------------------------------------------------------------
# Make a bar plot with error bars to eyeball levels of methylation across features
## -------------------------------------------------------------------------

#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 
## -------------------------------------------------------------------------

# Get the mean and confidence intervals for each feature per sample
summary_all<-summarySE(all, measurevar = "weightedMeth", 
                       groupvars = c("feature","sample"))

# Make a nice graph
ggplot(summary_all, aes(x=feature, y=weightedMeth, fill=sample))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth-ci, ymax=weightedMeth+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Weighted Methylation Level")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("sample1","sample2","sample3"),
                    labels=c("Sample 1","Sample 2","Sample 3"),
                    values=c("#44AA99","#CC6677","#DDCC77"))
