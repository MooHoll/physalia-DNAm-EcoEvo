# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

## Differential DNA methylation workshop Day 4 // 23rd April, 2026
### Instructor: James Ord

This is the Readme file for the session on Thursday, where we will cover differential DNA methylation analysis in `R`. All the data you'll need for this will be stored in: 

https://drive.google.com/drive/folders/1YILPV5wAoAwAh4kkBHpVl0QkDSEPOsJ5

The R packages we need are:
```
library(methylKit) # install via BiocManager::install()
library(matrixStats)
library(ggplot2)
library(biomaRt) # install via BiocManager::install()
library(genomation) # install via BiocManager::install()
```
And optionally:
```
library(pheatmap) for heatmap
library(dplyr) # for manhattan-like plot
#NOTE: be careful using the dplyr package alongside methylKit as both have functions called 'select()'
library(clusterProfiler) # for GO term enrichment, if desired; install via BiocManager::install()
```
## Overview

In this workshop we will quantify DNAme and calculate differential DNA methylation from RRBS data from threespine stickleback (G. aculeatus). The RRBS derived from gill tissue derived from either marine or freshwater / lake populations. The data were originally published by Artemov et al (2017) (https://academic.oup.com/mbe/article/34/9/2203/3813257). Sticklebacks are an interesting system to study adaptive evolution and plasticity, as they have repeatedly colonised freshwater habitats from ancestral marine environments, resulting in distinct marine and lake morphs, the latter being adapted to lower salinity. Whether epigenetic variation plays a role in this adaptation has been a question of gret interest. We will use a subset of the data to facilitate two comparisons: (1) a 'population' comparison between marine vs freshwater fish each kept at their native salinities (marine fish kept in seawater vs freshwater fish kept in lake water), and (2) an experimental comparison between marine fish in native salinity vs marine fish kept in lake water.

We will start with methylation counts files (.cov.gz) output from bismark_methylation_extractor after aligning RRBS reads to the stickleback V5 genome assembly obtained from the ENSEMBL database (https://www.ensembl.org/Gasterosteus_aculeatus/Info/Index). All sites are reference CpG context. Sites were not SNP-filtered. We will analyse these files using the R package methylKit, a popular package for differential methylation analyses.

An additional script is provided to process the same data using a different package, DSS, although methylKit has the option to implement the same statistical model for DMC / DMR calling as used by DSS (Bayesian beta binomial model).

## Workflow

### Reading in files / initial processing
Download the folder 'cov_files' from the google drive into your working directory.

If you check the top of just one of these files, e.g.
```
head(read.table(gzfile("cov_files/SRR3632630.CpG_report.merged_CpG_evidence.cov.gz")))
```
...you'll recognise the format of the cov files produced on Day 2 for the Arabidopsis data, with columns for chromosome, start, end, % methylation, number of Cs and number of Ts.
```
  V1   V2   V3       V4 V5 V6
1  I 2564 2565 95.23810 40  2
2  I 2691 2692 76.92308 10  3
3  I 2719 2720 92.30769 12  1
4  I 2725 2726 84.61539 11  2
5  I 2728 2729 92.30769 12  1
6  I 2733 2734 84.61539 11  2
```
We can read all of the cov files into R with methylKit's methRead() function. Note here that we specify that the counts files derived from bismark, and that we set a minimum coverage of 5 at each site (which is set here mostly to prune the data for performance purposes).
As for the sample IDs, the first three samples are 'MM' (marine fish in salt water), the next three are 'MF' (marine fish in freshwater) and the last three are 'FF' (freshwater fish in freshwater). By default, the samples are entered in alphabetical order of filename, and we specify the sample IDs and corresponding treatment codes in the order in which the files are read in.
```
#set path to methylation files
meth_path<-"cov_files/"
# make a list of the files
file.list<-as.list(list.files(meth_path,pattern = "\\.cov.gz$",full.names = TRUE))

# Generate the initial methylKit object with methRead
myobj<-methRead(file.list,
                sample.id=as.list(c("MM1","MM3","MM4",
                                  "MF3","MF4","MF5",
                                  "FF2","FF3","FF4")),
                assembly="stickleback_v5",
                treatment=c(0,0,0,1,1,1,2,2,2),
                context="CpG",
                mincov = 5, # this can reasonably be set lower for DMR analyses; set at 5 here mainly for performance purposes
                pipeline = "bismarkCoverage"
)
```
OPTIONAL: If you want to reduce the computation time for testing purposes, you can subset the data for just a few chromosomes. This speeds things up especially at the tileMethylCounts() step.
```
my_chroms <- c("I", "II", "III", "IV")
subset_list <- lapply(myobj, function(x) x[x$chr %in% my_chroms,])
subset_list <- as(subset_list, "methylRawList")
subset_list@treatment<-c(0,0,0,1,1,1,2,2,2) # we need to add the treatment variable back in here
myobj<-subset_list
```

methylKit allows the identification either of individual differentially methylated cytosines (DMCs) Or differentially methylated regions (DMRs); the pipeline is almost identical, the difference being that we summarise the counts into windows for DMRs. This means that for DMRs, it is not strictly necessary to set a minimum coverage threshold because counts at low coverage sites will be added to a total count within the window. However it is a good idea to filter out sites with very high coverage (99.9th percentile) as such abnormally high coverage is often the result of repetitive elements. The filterByCoverage() allows us to apply this filter:
```
# Filtering options for DMC identification
filtered.myobj=filterByCoverage(myobj,lo.count=NULL,lo.perc=NULL,
                                 hi.count=NULL,hi.perc=99.9)
```
Next, tileMethylCounts() summarises the counts into windows - this is a key step for DMR analysis but would be omitted for DMC analysis.

NOTE: It can take a while for methylKit to perform this step (up to 10 mins on my mac). You may wish to consider running this on a subset of the data first, to make sure it works (see subsetting step). If you are using mac or linux, you can increase 'mc.cores' to a slightly higher number (e.g. 6) to speed things up, depending on how many cores are available on your system.
```
tiles<-tileMethylCounts(filtered.myobj,win.size=1000,step.size=1000,cov.bases = 10, mc.cores=1)
```
ALSO NOTE: if you have pre-defined regions of interest (e.g. promoters, CpG islands), you may consider the regionCounts() method which will obtain total counts of methylated and unmethylated cytosines within each pre-defined feature of interest (see methylKit manual for this).

Some initial plotting of methylation and coverage distributions of individual samples:
```
# plot distributions of one sample
par(mfrow=c(2,2))
# for individual sites:
getMethylationStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE)
# and for the tiling windows:
getMethylationStats(tiles[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(tiles[[1]],plot=TRUE,both.strands=FALSE)

# plot distributions of all samples (for tiling windows only):
par(mfrow=c(3,3))
for (i in seq(1,9)){
  getMethylationStats(tiles[[i]],plot=TRUE,both.strands=FALSE)}
par(mfrow=c(3,3))
for (i in seq(1,9)){
  getCoverageStats(tiles[[i]],plot=TRUE,both.strands=FALSE)}

```
unite() will drop windows that are not covered in all samples
```
# merge samples
meth<-unite(tiles)
# then, optionally - nullify objects we no longer need, to save memoty
myobj<-NULL; filtered.myobj<-NULL; tiles<-NULL
```
sometimes we may want to remove chromosomes (e.g. sex chromosomes, organelles)
```
# see which chromosomes are represented
levels(as.factor(meth$chr))
# We will leave in the Y chromosome and as all the individuals are reportedly male, but if we wanted to
#meth<-meth[meth$chr!="Y", ]
# chromosome XIX also harbours a major sex determining region, which some people prefer to exclude.
# We may however want to remove unplaced contigs: these often contain highly repetitive sequences where many reads may not have properly mapped.
`%not_like%`<-purrr::negate(`%like%`)
meth<-meth[meth$chr%not_like%"JACDQR", ]
```
### Overview of methylation levels and variance
We saw above that one can plot the distributions of methylation levels for individual samples in the non-united object. We can also visualise means and variances using methylation levels extracted from the united object.
```
# get the methylation levels for each sample and window
pm<-percMethylation(meth)
head(pm)

# get the mean methylation leves of all sites or windows
methmeans<-rowMeans(pm)

# what about a specific group? Here the 'FF' samples are in columns 7,8, and 9 of the matrix
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
What do you notice about the spread of the samples across the two principal components? Are there any interesting overall patterns? Do the groups seem to have equal variance?

In case you want to plot the data yourself, e.g. with ggplot, you can extract the eigenvector and eigenvalues:
```
PCA.obj<-PCASamples(meth,obj.return = T)
# get the PC values for each sample
PCA.obj$x
# get the % variance explained for the first two PCs
PCA.obj$sdev[c(1,2)]/sum(PCA.obj$sdev)
```
### Calculating differential methylation
methylKit differential methylation analyses are limited to two-group comparisons, although covariates can be included (see the methylKit vignette). For models with multiple groups or factorial designs, check out DSS and / or edgeR.

Let's perform a differential methylation comparison between marine and freshwater fish. We first need to create a new methylKit object with only these six samples.
```
meth_pops <- reorganize(meth,sample.ids=c("MM1","MM3","MM4","FF2","FF3","FF4"),
                           treatment=c(0,0,0,1,1,1))
```
Note that here we specify the M group as control (0) and the freshwater group  as 'test' (1), so hypermethylation will correspond with increased methylation in the freshwater population and hypomethylation will correspond with drecreased methylation in freshwater.
How would you run the above command for the other comparison (MM vs MF)?

The default workflow employs the calculateDiffMeth() command to perform the differential methylation calculations using the logistic regression method. However, in this case we will use the Bayesian beta binomial method (originally implemented in the DSS package, and implemented in methylKit as the calculateDiffMethDSS() command), as it is able to handle unequal variance between groups. This will produce the methylation difference and q-value (adjusted p-value) for each region (or site, if performing site-level analysis). 
```
Diff_pops<-calculateDiffMethDSS(meth_pops)
head(Diff_pops)
```
Diff_pops is now a methylKit object that contains information pertaining to differential methylation. We can then immediately make a plot to visualise the % of windows on each chromosome that are differentially methylated:
```
diffMethPerChr(Diff_pops,plot=T,qvalue.cutoff=0.05, meth.cutoff=25)
# or just the numbers...
diffMethPerChr(Diff_pops,plot=F,qvalue.cutoff=0.05, meth.cutoff=25)
```
We can further use methylKit functions to narrow down on the DMRs themselves.
```
# e.g., get all DMRs with criteria of 25% difference in methylation level at q < 0.01
Diff_pops_25p<-getMethylDiff(Diff_pops,difference=25,qvalue=0.05)
Diff_pops_25p
# If you wish, you can also select for either hyper or hypomethylated regions:
# get only hypermethylated (increased methylation in freshwater):
Diff_pops_25p.hyper<-getMethylDiff(Diff_pops,difference=25,qvalue=0.05,type="hyper")
Diff_pops_25p.hyper # how many do we have?
# get only hypomethylated (decreased methylation in freshwater):
Diff_pops_25p.hypo<-getMethylDiff(Diff_pops,difference=25,qvalue=0.05,type="hypo")
Diff_pops_25p.hypo # how many do we have?
```
You can also extract the results as a dataframe to work with more flexibly...
```
DMdata<-getData(Diff_pops)
head(DMdata)
# classify DMRs according to your own criteria.
DMdata$result<-ifelse(DMdata$meth.diff>25&DMdata$qvalue<0.05,"hyper",
                      ifelse(DMdata$meth.diff< -25&DMdata$qvalue<0.05,"hypo","non_DM"))
# see numbers of regions in each category
nrow(subset(DMdata,result=="hypo"))
nrow(subset(DMdata,result=="hyper"))
nrow(subset(DMdata,result=="non_DM"))
```
For example, this allows us to write out BED files: we can use these to see where are the DMRs on a genome browser (e.g. IGV).
```
write.table(data.frame(chr=subset(DMdata,result!="non_DM")$chr,
                       start=subset(DMdata,result!="non_DM")$start-1,
                       end=subset(DMdata,result!="non_DM")$end-1,
                       result=subset(DMdata,result!="non_DM")$result),
            file="DMRs.bed",row.names = F,col.names = F,quote = F,sep="\t")

# note that we subtract 1 from start end coordinates, so that they become '0-based' coordinates (as per BED format)
```

### More visualisation

A nice way to visualise the overall landscape of DMRs is to plot the % meth. diff. of all the windows across a chromosome. Below make a plot for one chromosome, colouring the points according to whether they are hypermethylated, hypomethylated, or not differentially methylated.

```
# make a plot of meth. differences on chromosome 1
ggplot(subset(DMdata,chr=="I"),aes(x=start,y=meth.diff))+
         geom_point(aes(fill=result),pch=21,alpha=0.5)+
  scale_fill_manual(values=c("red","blue","grey"))+
  geom_hline(yintercept=c(-25,25),linetype="dashed")
```
Do you notice anything striking about the overall patterns on one chromosome?

NOTE: We can also make something akin to a Manhattan plot, showing the DMR distributions across all chromosomes on the same panel. However the code is a bit long, so this is included as a 'bonus' at the end.

A basic heatmap can also be made. This works better with a smaller number of DMRs (e.g. with the biggest % meth differences or in larger windows).
The heatmap can be useful for further visualising the variation amongst samples. In this case we see that while the DMRs distinguish FF from MM samples, there is one FF sample that is relatively more similar to the MM samples.
```
# get the DMRs according to criteria and make an index based on the chromosome names and start pos.
DMdata_DM<-subset(DMdata,result!="non_DM" & abs(meth.diff)> 50 & qvalue<0.01)
DMdata_DM$index<-paste(DMdata_DM$chr,DMdata_DM$start,sep="_")

# get data frame of percentage meth per sample (derived from percMethylation(meth) earlier)
pm<-percMethylation(meth)
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

# make the heatmap
pheatmap(meth3mat,clustering_method = "complete", show_rownames=T)
```

### Intersecting with gene structure annotation
methylKit objects can be used with functions from the genomation package to incorporate gene structure annotation.

We have a BED12 file with exon annotations for each transcript:
```
head(read.table("Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.agat_edit.bed"))
```
readTranscriptFeatures() reads in the BED12 file and converts it to a GRangesList object
```
gene.obj=readTranscriptFeatures("Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.agat_edit.bed")
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
For each DMR, we now have the distance to the nearest TSS, and the ID of that transcript. Let's plot those distances:
```
hist(TSSdists$dist.to.feature)
```
We expect a fairly tight distribution here as this is RRBS data, so the DMRs should be mostly in CG-rich regions and in the vicinity of genes.

Note that the output from getAssociationWithTSS() does not contain the coordinates of the DMRs. To make something more intuitive, add the TSS information to our table of DMRs, so that for each DMR we will have its location in the genome, the distance to the nearest transcript TSS, and the ID of said transcript:
```
Diff_pops_25p$closest_feature<-TSSdists$feature.name
Diff_pops_25p$closest_feature_TSSdist<-TSSdists$dist.to.feature
Diff_pops_25p$closest_feature_strand<-TSSdists$feature.strand
```

We may be particularly interested in the transcripts that are nearest to the DMRs. To narrow down on these transcripts, let's subset those with TSS overlapping with possible promoter regions.

Here I arbitrarily define the promoter region as being from 2000bp upstream to 500bp downstream of the TSS:

```
DM_promoter_overlap<-subset(getData(Diff_pops_25p),closest_feature_TSSdist > -2000 & closest_feature_TSSdist < 500)
head(DM_promoter_overlap)
nrow(DM_promoter_overlap) # how many DMRs overlap with arbitrarily defined promoters?
```

### Adding more information from the ENSEMBL database

Given that we have ENSEMBL IDs for the transcripts in close vicinity to DMRs, we can query the ENSEMBL database to get more information about the corresponding genes, e.g. gene names and descriptions if they exist. This can also be used to extract information about homology to other species, GO terms, and more.

First get a nonredundant list of transcript IDs, as some 'promoters' may overlap with multiple DMRs, thus duplicate entries.
```
DM_promoters<-DM_promoter_overlap$closest_feature[!duplicated(DM_promoter_overlap$closest_feature)]
library(biomaRt)
```
Now we are ready to query the database! Below, we submit the list of transcript IDs as the query to BiomaRt, tellinf it that we want to pull down the gene id, gene name, and description for each one.

NOTE: the ENSEMBL servers may not always be responsive or operational. If you cannot connect, you can try a different mirror (see the help for useEnsembl() command), or just come back later and try again.

```
ensembl <- useEnsembl(biomart = "genes", dataset = "gaculeatus_gene_ensembl")
listAttributes(ensembl) # list possible attributes
# Query the database via BiomaRt to extract information
query<-getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id','external_gene_name', 'description'),
             values = DM_promoter_overlap$feature.name, 
             filters='ensembl_transcript_id',
             mart = ensembl)

head(query)
```

Note also: The annotations of the genome version used for this analysis are from release 115, the most recent at the time of writing. Ideally, you will always be working with the most up to date annotation.

Finally, let's merge the results from the query with our 'DM_promoter_overlap' dataframe. Then for each DMR overlapping putative promoters we will have the location, DMR statistics, closest TSS transcript and the name / description of the gene if available.

```
colnames(query)[1]<-"closest_feature"
DM_promoter_overlap<-merge(DM_promoter_overlap,query,by="closest_feature")
```

We can also write this out as a table, e.g. for browsing in excel:
```
write.table(DM_promoter_overlap,file="DMRs_promoter_overlap.txt",sep="\t",col.names = T,row.names = F,quote=F)
```

## Exercises:
* Perform the DMR analysis for the environmental comparison (MM vs MF)
* Do any transcripts have differentially methylated promoter regions in both comparisons?
* What is the relative distribution of hyper- and hypomethylated regions?
* Is there a difference in the distribution of covered features (exons, promoters etc.) between hypo- and hyper-DMRs?
* Do these two classes of DMRs appear to be related to different functions?

### BONUS CODE: Manhattan-like plot of differential methylation

This code uses the 'DMdata' dataframe generated just after the differential methylation calculation step.

```
# modified from code by Daniel Roelfs:
# https://danielroelfs.com/posts/how-i-create-manhattan-plots-using-ggplot/

# add the cumulative coordinates so that we can display points from all chromosomes sequentially
data_cumul <- DMdata |>
  group_by(chr) |>
  summarise(max_bp = max(start)) |>
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
  select(chr, bp_add)
DMdata <- DMdata |>
  inner_join(data_cumul, by = "chr") |>
  mutate(bp_cumul = start + bp_add)

# set the points for the x axis (chromosome centres)
axis_set <- DMdata |>
  group_by(chr) |>
  summarize(center = mean(bp_cumul))

# code the plot
manhplot <- ggplot(subset(DMdata,result=="non_DM"), aes(
  x = bp_cumul, y = meth.diff,
  color = as.factor(chr))) +
  geom_hline(
    yintercept = c(-25,25), color = "grey40",
    linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_x_continuous(
    label = axis_set$chr,
    breaks = axis_set$center) +
  scale_color_manual(values = rep(c("black", "darkgrey"),
    unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = NULL,y = "% methylation difference") +
  theme_minimal() +
  theme(legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))+
  geom_point(data=subset(DMdata,result=="hyper"),pch=21,mapping=aes(x = bp_cumul, y = meth.diff),fill="red")+
  geom_point(data=subset(DMdata,result=="hypo"),pch=21,mapping=aes(x = bp_cumul, y = meth.diff),fill="blue")
```

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
enricher(gene = DM_promoters,TERM2GENE = term2gene, TERM2NAME = term2name,universe = all_covered_transcripts,pAdjustMethod = "fdr", pvalueCutoff = 0.1)
```
Appendix: Obtaining BED12 file of stickleback transcripts
The stickleback V5 gff3 annotation was obtained (https://ftp.ensembl.org/pub/release-115/gff3/gasterosteus_aculeatus/) and converted it to BED12 format (for use with genomation) with the AGAT toolkit:
```
wget https://ftp.ensembl.org/pub/release-115/gff3/gasterosteus_aculeatus/Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.gff3.gz
agat_convert_sp_gff2bed.pl -gff Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.gff3.gz -o Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.agat.bed
sed -e 's/transcript://g' Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.agat.bed > Gasterosteus_aculeatus.GAculeatus_UGA_version5.115.agat_edit.bed
```
