# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

https://www.physalia-courses.org/

### Hollie Marshall, James Ord, Alexandre de Mendoza

#### 10th &ndash; 14th February, 2025

## Overview

Epigenetic studies on DNA methylation are now possible to do at the nucleotide level genome-wide and also in non-model organisms. In this course, we will introduce the different available approaches for obtaining and analysing DNA methylation data, including bisulfite sequencing (BS-seq, EM-seq) with Ilumina and long reads with Oxford Nanopore. We will also touch on PacBio long-read analysis. We will cover all necessary steps to obtain methylation information from high throughput data to statistical analyses used to identify differentially methylated sites and regions. The data will be interpreted in terms of their biological importance in the field of ecology and evolution.

## Target audience and assumed background

This course is aimed at researchers and technical workers who are or who will be generating and/or analysing DNA methylation from whole genome bisulfite sequencing data from IIlumina, PacBio data or Nanopore data. Examples demonstrated in this course will involve primarily non-model organisms with a draft reference genome available and examples of applications of this data type for different purposes will be covered.

## Program

### Monday: DNA methylation analysis with short reads (Illumina WGBS and EM-Seq) – 2-7pm Berlin time

1. Presentation: Origin and evolution of DNA methylation in eukaryotes: from parasite defence to gene regulation.
2. Presentation: Introduction to short-read format and processing (WGBS/EM-Seq).
3. Lab 1: Short-read processing:
* QC, filtering, trimming
* Genome preparation
* Alignment, deduplication

### Tuesday: DNA methylation analysis with short reads (Illumina WGBS and EM-Seq) – 2-7pm Berlin time

1. Presentation: Introduction to DNAm in plants.
2. Presentation: Population epigenomics: patterns and drivers of DNAm variation in natural populations.
3. Lab 2: Short-read processing:
* Conversion efficiencies (lambda/pUC19)
* Methylation count files
* Accounting for SNPs

### Wednesday: DNA methylation analysis with long reads (Oxford Nanopore) – 2-7pm Berlin time

1. Presentation: Nuances of DNA methylation data: sex, tissue, age, imprinting and allelic-variation.
2. Presentation: Introduction to long-read format and processing (Nanopore/PacBio).
3. Lab 3: Nanopore long-read processing:
* Fast5/pod5 generation, Pod5 tools
* Dorado/Guppy, and models
* Modbam2bed/Modkit
* Example: methylation haplophasing

### Thursday: Identification of differentially methylated sites and regions and integration with SNP data – 2-7pm Berlin time

1. Presentation: DNAm and genetic variation: biological insights and analytical hazards.
2. Presentation: Introduction to differential DNA methylation.
3. Lab 4: Differential DNA methylation:
* Testing for differentially methylated sites and regions between populations and conditions
* Visualising differential methylation
* Accounting for SNPs

### Friday: Integration of DNA methylation and gene expression and open session on your own data / experiments – 2-7pm Berlin time
1. Lab 5: Integration of DNA methylation and gene expression:
* Quantifying genes/exon/intron levels
* Correlational plots
* Overlap of differential methylation and differential expression
2. Open session: Group questions / discussion:
* Follow up on any course material
* Student questions on their own data / experimental plans

## Pre-course preparation

### R

Please insure you have version 4.4.x of R installed. Note that R and RStudio are two different things:
it is not sufficient to just update RStudio, you also need to update R by
installing new versions as they are released.

To download R go to the [CRAN Download](https://cran.r-project.org/) page and 
follow the links to download R for your operating system:

* [Windows](https://cran.r-project.org/bin/windows/)
* [MacOS X](https://cran.r-project.org/bin/macosx/)
* [Linux](https://cran.r-project.org/bin/linux/)

To check what version of R you have installed, you can run

```r
version
```

in R and look at the `version.string` entry (or the `major` and `minor`
entries).

### RStudio

We will use RStudio as the environment for interacting with R. You are free to
use your preferred IDE or R interface if you wish. To install RStudio, go to
<https://posit.co/download/rstudio-desktop/> and follow the instructions.

