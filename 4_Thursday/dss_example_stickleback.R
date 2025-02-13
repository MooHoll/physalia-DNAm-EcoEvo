#BiocManager::install("DSS")
library(DSS)
require(bsseq)

####################
# Core DSS worflow #
####################

# Read in all files
file.list = list.files("./stickleback_cov_files/",full.names = T)

# The data needs to be in a slightly specific format
# Unlike the default .cov format read by bismark, we have columns for chr, pos, coverage (N) and methyl count (X)
# The following function reads a file and reformats accordingly
read_and_reformat <- function(x){
  y<-read.table(gzfile(x))
  y$N<-y$V5+y$V6
  y<-y[c(1,2,7,5)]
  colnames(y)<-c("chr","pos","N","X")
  #y<-subset(y,chr=="chrI") # un-comment to select just one chromosome for testing purposes
  return(y)
  }

# example with one file
head(read_and_reformat("./stickleback_cov_files_130225/FF1.merged_CpG_evidence.SNPfilt.mincov5.cov.gz"))

# now use lapply to apply the function to all the files
samples <- lapply(file.list, read_and_reformat)

# the reformated sample files are now stored in 'samples'
# now we process them with DSS as follows. Note the order of sample names.
BSobj = makeBSseqData(samples,
                      c("FF1","FF2","FF3",
                        "MF1","MF2","MF3",
                        "MM1","MM2","MM3"))
BSobj

# Note that no filtering is performed on coverage: DSS attempts to use all of the information even from lowly-covered sites.

# here we test for differential methylation between two groups
# this runs the smoothing algorithm and estimates dispersion at the CpG site level
# this usually takes some time
dmlTest.sm1 = DMLtest(BSobj, group1=c("MM1","MM2","MM3"), group2=c("FF1","FF2","FF3"), 
                     smoothing=TRUE)

# call DMCs (or DMLs as they are called in this package) from the dmlTest result
dmls1 = callDML(dmlTest.sm1, delta=0.25, p.threshold=0.01)
head(dmls1)

# call DMRs using two thresholds
dmrs1 = callDMR(dmlTest.sm1, delta=0.25, p.threshold=0.05,minCG = 5)
dmrs2 = callDMR(dmlTest.sm1, delta=0.25, p.threshold=0.05,minCG = 10)
head(dmrs1)

# Plot one of the differentially methylated regions
showOneDMR(dmrs1[2,], BSobj)

####################################
# Overlapping with gene annotation #
####################################

library(GenomicRanges)

# read in the transcript annotation BED file
transcript_annot<-read.table("stickleback_v5_ensembl_genes.bed")
head(transcript_annot)

# Here we convert the transcript annotation BED file to a GRanges object, where "seqnames" is the chromosome name column
a <- GRanges(seqnames = transcript_annot$V1, ranges = IRanges(start=transcript_annot$V2-999, end=transcript_annot$V3+1001,mcols=promoters$V4))
b <- GRanges(seqnames = dmrs1$chr, ranges = IRanges(start=dmrs1$start, end=dmrs1$end,mcols=dmrs1$diff.Methy))
overlap_promoter_transcripts<-subsetByOverlaps(a,b)

# see how many promoters overlap with DMRs:
length(overlap_promoter_transcripts$mcols)

# but which genes correspond with which DMRs?
fo<-findOverlaps(a,b)
fo
# this shows us that, e.g. the 23rd, 24th, and 25th DMRs all overlap with transcript 604, which is:
transcript_annot[604,]
# show the DMR:
showOneDMR(dmrs1[23,], BSobj)

# We could now use the ENSEMBL biomaRt package to extract annotations for these transcripts, as we did for the methylKit example
