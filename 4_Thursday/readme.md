# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

## Differential DNA methylation workshop Day 4 // 13th February, 2025
### Instructor: James Ord

This is the Readme file for the session on Thursday, where we will cover differential DNA methylation analysis in `R`. All the data you'll need for this will be stored in: 

`~/Share/4_Thursday/`

## Overview

In this workshop we will quantify DNAme and calculate differential DNA methylation from RRBS data from threespine stickleback (G. aculeatus). The RRBS derived from gill tissue derived from either marine or freshwater / lake populations. The data were originally published by Artemov et al (2017) (https://academic.oup.com/mbe/article/34/9/2203/3813257) and subsequently re-analysed by Ord et al (2023) (https://academic.oup.com/mbe/article/40/4/msad068/7083720). Sticklebacks are an interesting system to study adaptive evolution and plasticity, as they have repeatedly colonised freshwater habitats from ancestral marine environments, resulting in distinct marine and lake morphs, the latter being adapted to lower salinity. Whether epigenetic variation plays a role in this adaptation has been a question of gret interest. We will use a subset of the data to facilitate two comparisons: (1) a 'population' comparison between marine vs freshwater fish each kept at their native salinities (marine fish kept in seawater vs freshwater fish kept in lake water), and (2) an experimental comparison between marine fish in native salinity vs marine fish kept in lake water.

We will start with methylation counts files (.cov.gz) output from bismark_methylation_extractor after aligning RRBS reads to the stickleback V5 genome assembly (https://stickleback.genetics.uga.edu/downloadData/). All sites are reference CpG context. The counts files were SNP-filtered using a combination of SNP-callers performed on the RRBS data themselves (BS-SNPer, Biscuit, and CGmaptools). The non-SNP-filtered counts files are also provided, however. We will analyse these files using the R package methylKit, a popular package for differential methylation analyses.

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
# optional - nullify objects we no longer need
myobj<-NULL; tiles<-NULL
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
Diff_pops_25p # reveals 1,639 DMRs at these criteria
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
### Intersecting with functional annotation
methylKit objects can be used with functions from the genomation package to incorporate functional annotation.
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
* What happens if you use non-SNP-filtered counts files? does this affect the number of DMRs found?

Appendix: Obtaining BED12 file of stickleback transcripts
The stickleback V5 gff3 annotation was obtained (https://stickleback.genetics.uga.edu/downloadData/) and converted it to BED12 format (for use with genomation) with the AGAT toolkit:
```
agat_convert_sp_gff2bed.pl -gff stickleback_v5_ensembl_genes.gff3 -o stickleback_v5_ensembl_genes.bed
```
