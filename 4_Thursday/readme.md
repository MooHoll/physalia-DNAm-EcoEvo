# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

## Differential DNA methylation workshop Day 4 // 13th February, 2025
### Instructor: James Ord

This is the Readme file for the session on Thursday, where we will cover differential DNA methylation analysis in `R`. All the data you'll need for this will be stored in: 

https://drive.google.com/drive/folders/1cGRQGoylbbWUzECj_ATLJIyHoBPCwNOn?usp=sharing

The R packages we need are:
```
library(methylKit)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(genomation)
library(GenomicRanges)
```
And optionally:
```
library(clusterProfiler) # for GO term enrichment, if desired
library(DSS) # for alternative DMC/DMR calling with DSS
require(bsseq)
```
## Overview

In this workshop we will quantify DNAme and calculate differential DNA methylation from RRBS data from threespine stickleback (G. aculeatus). The RRBS derived from gill tissue derived from either marine or freshwater / lake populations. The data were originally published by Artemov et al (2017) (https://academic.oup.com/mbe/article/34/9/2203/3813257) and subsequently re-analysed by Ord et al (2023) (https://academic.oup.com/mbe/article/40/4/msad068/7083720). Sticklebacks are an interesting system to study adaptive evolution and plasticity, as they have repeatedly colonised freshwater habitats from ancestral marine environments, resulting in distinct marine and lake morphs, the latter being adapted to lower salinity. Whether epigenetic variation plays a role in this adaptation has been a question of gret interest. We will use a subset of the data to facilitate two comparisons: (1) a 'population' comparison between marine vs freshwater fish each kept at their native salinities (marine fish kept in seawater vs freshwater fish kept in lake water), and (2) an experimental comparison between marine fish in native salinity vs marine fish kept in lake water.

We will start with methylation counts files (.cov.gz) output from bismark_methylation_extractor after aligning RRBS reads to the stickleback V5 genome assembly (https://stickleback.genetics.uga.edu/downloadData/). All sites are reference CpG context. The counts files were SNP-filtered using a combination of SNP-callers performed on the RRBS data themselves (BS-SNPer, Biscuit, and CGmaptools). The non-SNP-filtered counts files are also provided, however. We will analyse these files using the R package methylKit, a popular package for differential methylation analyses.

An additional script is provided to process the same data using a different package, DSS.

## Workflow

### Reading in files / initial processing
Download the folder 'Stickleback_covfiles_SNPfilt' from the server into your working directory, either using a file transfer client like filezilla or using the sftp command on a local terminal.
Once they are accessible on your local computer, we can read them into R with methylKit's methRead() function. Note here that we specify that the counts files derived from bismark, and that we set a minimum coverage of 5 at each site.
```
#set path to methylation files
meth_path<-"Stickleback_covfiles_SNPfilt/"
# make a list of the files
file.list<-as.list(list.files(meth_path,pattern = "\\.cov.gz$",full.names = TRUE))

# Generate the initial methylKit object with methRead
myobj<-methRead(file.list,
                sample.id=as.list(c("FF1","FF2","FF3",
                                  "MF1","MF2","MF3",
                                  "MM1","MM2","MM3")),
                assembly="stickleback_v5",
                treatment=c(1,1,1,2,2,2,0,0,0),
                context="CpG",
                mincov = 5,
                pipeline = "bismarkCoverage"
)
```
methylKit allows the identification either of individual differentially methylated cytosines (DMCs) Or differentially methylated regions (DMRs); the pipeline is almost identical, the difference being that we summarise the counts into windows for DMRs. Additional filtering may be desired for site-level analysis, e.g.:
```
# Filtering options for DMC identification
filtered.myobj=filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                                 hi.count=NULL,hi.perc=99.9)
```
tileMethylCounts() summarises the counts into windows
```
tiles<-tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 5)
```
Some initial plotting of methylation distributions of individual samples:
```
# plot distributions of one sample
getMethylationStats(tiles[[2]],plot=TRUE,both.strands=FALSE)

# plot distributions of all samples
par(mfrow=c(3,3))
for (i in seq(1,9)){
  getMethylationStats(tiles[[i]],plot=TRUE,both.strands=FALSE)}

```
unite() will drop windows that are not covered in all samples
```
# merge samples
meth<-unite(tiles)
```
sometimes we may want to remove chromosomes (e.g. sex chromosomes, organelles)
```
# see which chromosomes are represented
levels(as.factor(meth$chr))
# Remove the Mitochondrion; We will leave in the Y chromosome and as all the individuals are reportedly male
meth<-meth[meth$chr!="chrM", ]
```
### Overview of methylation levels and variance
We saw above that one can plot the distributions of methylation levels for individual samples in the non-united object. We can also visualise means and variances using methylation levels extracted from the united object.
```
# get the methylation levels for each sample and window
pm<-percMethylation(meth)
head(pm)

# get the mean methylation leves of all sites or windows
methmeans<-rowMeans(pm)

# what about a specific group? Here the 'MM' samples are in columns 7,8, and 9 of the matrix
methmeans_MM<-rowMeans(pm[,c(7,8,9)])
methSDs_MM<-rowSds(pm[,c(7,8,9)])
# plot the mean and SD
par(mfrow=c(2,1))
hist(methmeans_MM)
hist(methSDs_MM)
```
How do means and SDs compare between groups?

It may be pertinent to exclude sites that have no variation, as they are not informative. We can use our array of mean methylation levels to filter the methylation object to remove those windows with no variation.
```
meth<-meth[methmeans>0&methmeans<100,]
```
### Futher exploratory analysis: PCA
PCA is a common way of inspecting the extent of variation within and between groups. In methylKit we can generate and plot this with one command.
```
par(mfrow=c(1,1))
PCASamples(meth)
```
In case you want to plot the data yourself, e.g. with ggplot:
```
PCA.obj<-PCASamples(meth,obj.return = T)
# get the PC values for each sample
PCA.obj$x
# get the % variance explained for the first two PCs
PCA.obj$sdev[c(1,2)]/sum(PCA.obj$sdev)
```
### Differential methylation
methylKit differential methylation analyses are limited to two-group comparisons, although covariates can be included (see the methylKit vignette). For models with multiple groups or factorial designs, check out DSS and / or edgeR.

Let's perform a differential methylation comparison between marine and freshwater fish. We first need to create a new methylKit object with only these six samples.
```
meth_pops <- reorganize(meth,sample.ids=c("MM1","MM2","MM3","FF1","FF2","FF3"),
                           treatment=c(0,0,0,1,1,1))
```
Next, calculateDiffMeth() will actually perform the logistic-regression based calculations of differential methylation. This will produce the methylation difference and q-value (equivalent to an adjusted p-value) for each region (or site, if performing site-level analysis):
```
Diff_pops<-calculateDiffMeth(meth_pops)
head(Diff_pops)

# A plot to visualise the % of windows on each chromosome that are differentially methylated
diffMethPerChr(Diff_pops,plot=T,qvalue.cutoff=0.05, meth.cutoff=25)
# or just the numbers...
diffMethPerChr(Diff_pops,plot=F,qvalue.cutoff=0.05, meth.cutoff=25)
```
Diff_pops is now a methylKit object that contains information pertaining to differential methylation. We can further use methylKit functions to narrow down on the DMRs themselves.
```
# e.g., get all DMRs with criteria of 25% difference in methylation level at q < 0.01
Diff_pops_25p<-getMethylDiff(Diff_pops,difference=25,qvalue=0.05)
Diff_pops_25p
# get only hypermethylated (increased methylation in freshwater):
Diff_pops_25p.hyper<-getMethylDiff(Diff_pops,difference=25,qvalue=0.05,type="hyper")
```
Alternatively, you can extract the results as a dataframe to explore more flexibily...
```
DMdata<-getData(Diff_pops)
head(DMdata)
# classify DMRs according to your own criteria.
DMdata$result<-ifelse(DMdata$meth.diff>25&DMdata$qvalue<0.05,"hyper",
                      ifelse(DMdata$meth.diff< -25&DMdata$qvalue<0.05,"hypo","non_DM"))
# see numbers of regions
nrow(subset(DMdata,result=="hypo"))
nrow(subset(DMdata,result=="hyper"))
nrow(subset(DMdata,result=="non_DM"))

# make a plot of meth. differences on chromosome 1
ggplot(subset(DMdata,chr=="chrI"),aes(x=start,y=meth.diff))+
         geom_point(aes(fill=result),pch=21,alpha=0.5)+
  scale_fill_manual(values=c("red","blue","grey"))
```
A basic heatmap can also be made. This works best if you have a smaller number of DMRs (e.g. with the biggest % meth differences or in larger windows.
```
# HEATMAP OF SOME DMRs
library(pheatmap)
# get the DMRs according to criteria and make an index based on the chromosome names and start pos.
DMdata_DM<-subset(DMdata,result!="non_DM" & abs(meth.diff)> 25 & qvalue<0.01)
DMdata_DM$index<-paste(DMdata_DM$chr,DMdata_DM$start,sep="_")

# get data frame of percentage meth per sample (derived from percMethylation(meth) earlier)
meth2<-as.data.frame(pm)
# add an index from the object 'DMdata' (derived from getData(Diff_pops) earlier)
# --> the rows should be in the same order in DMdata as they are in the & meth dataframe
meth2$index<-paste(DMdata$chr,DMdata$start,sep="_")
rownames(meth2)<-meth2$index # replace rownames with the index
# (this is useful if you want to show the chrom_pos on the heatmap,
# but only works well if you have a small number of DMRs)

# now we can subset the % meth dataframe for only the indexes corresponding with DMRs
meth3<-subset(meth2,index %in% DMdata_DM$index)[-10] # last part removes the index column
# convert to matrix
meth3mat<-as.matrix(meth3)

pheatmap(meth3mat,clustering_method = "complete", show_rownames=F)
```

### Intersecting with gene structure annotation
methylKit objects can be used with functions from the genomation package to incorporate gene structure annotation.
```
library("genomation")
```
We have a BED12 file with exon annotations for each transcript:
```
head(read.table("stickleback_v5_ensembl_genes.bed"))
```
readTranscriptFeatures() reads in the BED12 file and converts it to a GRangesList object
```
gene.obj=readTranscriptFeatures("stickleback_v5_ensembl_genes.bed")
gene.obj
```
Intersect the DMRs with the annotation. There may be some warnings about chromosomes not present in either DMR set or annotation
```
Diff_pops_25p_Ann<-annotateWithGeneParts(as(Diff_pops_25p,"GRanges"),gene.obj)
# See a summary of gene / promoter overlap:
getTargetAnnotationStats(Diff_pops_25p_Ann,percentage=TRUE)
# Same but as a plot:
plotTargetAnnotation(Diff_pops_25p_Ann,main="DMRs: marine vs freshwater")
```
Note that arbitrarily, promoters are defined as +/- 1kb from the transcription start site (TSS), but this can be changed in the options for readTranscriptFeatures()

To find the specific genes that are associated with our DMRs, Genomation has a neat function to link each feature (e.g. a DMR) with the TSS of the nearest gene:
```
TSSdists<-getAssociationWithTSS(Diff_pops_25p_Ann)
head(TSSdists)
```
For each DMR, we now have the distance to the nearest TSS, and the ID of that gene. Let's plot those distances:
```
hist(getAssociationWithTSS(Diff_pops_25p_Ann)$dist.to.feature)
```
We may be particularly interested in the transcripts that are nearest to the DMRs. To narrow down on these transcripts, let's subset those with TSS within 1kb of the DMRs (i.e. overlapping with an arbitrarily defined promoter region).
```
DM_promoters<-subset(TSSdists,dist.to.feature > -1000 & dist.to.feature < 1000)
head(DM_promoters)
nrow(DM_promoters) # putative promoters of 240 transcripts overlap DMRs
```

### Adding more information from the ENSEMBL database
Now we have a list of transcripts that are interesting to us because of DMRs overlapping the promoters. We can fetch more information, e.g. gene names and descriptions, from the ENSEMBL database...
```
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "gaculeatus_gene_ensembl",host="https://sep2019.archive.ensembl.org")
```
The annotations of the genome version used for this analysis are from release 95. Here we therefore use a 'legacy' version of the ENSEMBL database, release 99. Ideally, you will always be working with the most up to date annotation.
```
listAttributes(ensembl) # list possible attributes

# Query the database via BiomaRt to extract information
query<-getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id','external_gene_name', 'description'),
             values = DM_promoters$feature.name, 
             filters='ensembl_transcript_id',
             mart = ensembl)

head(query)
```
Note that there are 4 transcript entries missing. I suspect a small number were deprecated between ENSEMBL versions 95 and 99.

## Exercises:
* Perform the DMR analysis for the environmental comparison (MM vs MF)
* Are any promoter regions differentially methylated in both comparisons?
* What is the relative distribution of hyper- and hypomethylated regions?
* Is there a difference in the distribution of covered features (exons, promoters etc.) between hypo- and hyper-DMRs?
* Do these two classes of DMRs appear to be related to different functions?
* How do parameters such as site coverage and number of CpGs in the windows influence DMR detection?

### OPTIONAL: GO term enrichment analysis

Enrichment analysis tests for whether certain functional categories (e.g. GO terms) are enriched amongst gene lists (e.g. differentially expressed or differentially methylated), compared to a background. The results can depend on many factors including the completeness / quality of functional annotation and the choice of background set. Ideally the background should be genes (or transcripts in this case) that are represented in the analysis (i.e. covered in the DMR analysis). Results of these can give some vague indication of what processes might be affected, but they should be taken with a pinch of salt (particularly as the functional annotation usually comes from homology to model organisms, and is usually not very extensive in nonmodels). I include this as a demonstration of how it is done using information from ENSEMBL annotation and with the clusterProfiler package, but this should not be considered an essential step.

```
library(clusterProfiler)
# For our background, we want the transcript IDs of all promoters covered by the differential methylation analysis
Diff_pops_Ann<-annotateWithGeneParts(as(Diff_pops,"GRanges"),gene.obj)
# see how many transcript promoters are covered
all_covered_transcripts<-levels(as.factor(Diff_pops_Ann@dist.to.TSS$feature.name))
length(all_covered_transcripts)

query_goterms<-getBM(attributes = c('ensembl_transcript_id','go_id','name_1006'),
                     values = all_covered_transcripts, 
                     filters='ensembl_transcript_id',
                     mart = ensembl)
head(query_goterms)

# subset for biological processes
query_goterms_bp<-subset(query_goterms,namespace_1003=="biological_process")

# set up the maps linking (1) go terms with go term names and (2) go terms with transcript IDs
term2name<-query_goterms_bp[c(2,3)]
term2gene<-query_goterms_bp[c(2,1)]
# the universe, a.k.a. background -> should be all the transcripts covered that also have some GO annotation

# run the enricher() function
enricher(gene = DM_promoters$feature.name,TERM2GENE = term2gene, TERM2NAME = term2name,universe = all_covered_transcripts,pAdjustMethod = "fdr", pvalueCutoff = 0.1)

# How do the results differ between hypo- and hypermethylated genes?
# Is anything enriched amongst promoters that are differentially methylated in both experimental and population comparisons?
```
Appendix: Obtaining BED12 file of stickleback transcripts
The stickleback V5 gff3 annotation was obtained (https://stickleback.genetics.uga.edu/downloadData/) and converted it to BED12 format (for use with genomation) with the AGAT toolkit:
```
agat_convert_sp_gff2bed.pl -gff stickleback_v5_ensembl_genes.gff3 -o stickleback_v5_ensembl_genes.bed
```
