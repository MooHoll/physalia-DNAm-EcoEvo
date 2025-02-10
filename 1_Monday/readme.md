# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

# Short-reads workshop Day 1 // 10th February, 2025
### Instructor James Ord

This is the Readme file for the session on Monday, where we will cover the first half of the short-read processing. All the data you'll need for this will be stored in: 

`~/Share/1_Monday/`

## Overview

In today's and tomorrow's sessions we will go through a basic version of bisulfite sequencing data processing using well-established tools, namely the Bismark suite.
Today we will cover the initial stages of trimming, and alignment. We will use data from the popular model plant Arabidopsis thaliana (thale cress) which has a high quality reference genome and an extensive catalogue of genomic and epigenomic data from over a thousand naturally-occurring varieties (accessions). Here we'll go through an example with a single-end read file.

## 1. Obtain raw reads and subsample

Find the data for today's session in the folder:

`~/Share/1_Monday/raw_data_for_alignment`

You will find four folders, each corresponding with a different A. thaliana accession.

```
$ ls
Col-0  Dog-4  Etna-2  Kyoto
```
Two (Col-0 and Kyoto) contain single-end data, while the other two (Dog-4 and Etna-2) contain paired-end data. Let's look inside the Col-0 folder.
```
$ ls Col-0
ftp.txt  SRR771698.fastq.gz
```
The file `SRR771698.fastq.gz` is the raw BSseq (single-end) read file in FASTQ format. The 'SRR' denotes that this is a sequencing run that was submitted to the NCBI sequence read archive. These read files are usually also backed up on the European Nucleotide Archive (ENA). For example, the present file can be found here:
https://www.ebi.ac.uk/ena/browser/view/SRR771698
The file `ftp.txt` is the FTP (file transfer protocol) link for this file. Any such file available on ENA can be downloaded given such a link, using the `wget` command.

Let's take a look at the top of the file. The FASTQ format used for BSseq reads is a standard format for reads produced by most short read sequencing platforms.
Here we use the `zcat`and `head` command to see the first 10 lines of the file. `zcat` is a special case of `cat` (short for concatenate) which can be used on gzip-compressed (.gz) files. By itself, `zcat` would print the whole file, so we pipe the output (`|`) to the `head` command to see just the first ten lines:
```
$ zcat SRR771698.fastq.gz | head
@SRR771698.1 DARWIN_4066_FC6290TAAXX:5:1:1352:1030/1
TATACCGCAACCGGATTTTAAAGGCGTAAGAATTGTATTTTTGTTAGAAGATATAAAGTTAAAGATTTATATGGATTTTGGTTAT
+
HDHDHGDE<4BB<9?GGGGGGFGEDGE<GE?FFEF@@GGGGGGDDAGBBGDGBGD<2D??B??D8CECCA8EFE8CCACDC@>CC
@SRR771698.2 DARWIN_4066_FC6290TAAXX:5:1:1745:1027/1
NAAAGTTAAAGATTTATATGGATTTTGGTTATATTATGAAAGTTTTGAGAAGTAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGA
+
%+,*+3000.<:<<<@@@@@8888888888@@@@@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@SRR771698.3 DARWIN_4066_FC6290TAAXX:5:1:1789:1034/1
TTTTTTTGATAAAGGAAAAGGAAATATATATTTATTGGATTAATTAGTGATTAGTAGATTATTATGAAAAGGAAATAGAATTGAT
```
You can find more information about the FASTQ format here:
https://en.wikipedia.org/wiki/FASTQ_format

This file contains nearly 40 million short reads. To process all of these in this session would be impractical for learning purposes, so we will take a subsample of them. Subsampling is an extremely useful strategy in genomic data analyses as it allows us to test and troubleshoot our analysis workflows quickly and with minimal resources.

Subsampling can be performed quickly using the `reformat.sh` script from the bbtools suite:
```
conda activate short_reads # activate the conda environment
reformat.sh in=~/Share/1_Monday/raw_data_for_alignment/Col-0/SRR771698.fastq.gz out=Col0_subsample.fastq.gz samplerate=0.1
```
Note that with the `in=` argument we can specify the direct path to the file, so we don't need to make a copy of it. The last argument is the sampling rate which here we set to 0.1 to get 10% of the reads, but you can set it lower if you want to speed up the subsequent steps. The command will generate `Col0_subsample.fastq.gz` in your working directory.

Subsampling can also be done on paired read files, ensuring that for every read that is sampled, the other member of the pair is also sampled. See the help file (run `reformat.sh --help`) or the script `subsample_PE.sh` in this repository.

### Core Tasks:
* Inspect a fastq.gz file using zcat and head
* Obtain a subsample of a single-end fastq.gz file using `reformat.sh`

### Advanced Tasks:
* Perform subsampling on paired-end read files using `reformat.sh`
* Find and obtain data of another accession from ENA and subsample from it

## 2. Quality trimming and QC with FASTP

Trimming short reads for adapter sequences and low-quality bases is an important step prior to downstream analyses, as low quality bases and adapter sequences may result in false positive methylation calls. Several tools are available for quality trimming including `trimmomatic`, `cutadapt`, `bbduk.sh`, and `fastp`. Here we will use the `fastp` which generates HTML QC reports (similar to FASTQC) in addition to trimming the reads. Conveniently, also detects the adapter sequence automatically.

Check the `fastp` documentation here to get information about the various options:
https://github.com/OpenGene/fastp

Run `fastp` on a single-end file as follows:
```
fastp --thread 2 --in1 Col0_subsample.fastq.gz --out1 Col0_subsample.trimmed.fastq.gz \
    -l 35 -h Col0_subsample_SE.fastp.html \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20
```
Above takes a single end file (Col0_subsample.fastq.gz) and outputs a trimmed version. Here, bases are cut from the front (5') and tail (3') end of the read if the phred33 quality score is below 20, which translates to trimming the base if there is more than a 1% change of an incorrect base call. See more about phred33 scores here:
https://en.wikipedia.org/wiki/Phred_quality_score

You can adjust the quality parameters as you see fit. For the paired end version, see `fastp_PE.sh` in this repository.

In addition to producing the trimmed read file(s), fastp will also produce a report in HTML format. To view this you will need to download it to your local computer. This can be easily done in drag-and-drop fashion if you are using a linux client with FTP interface such as MobaXterm. Otherwise, you can open a new linux terminal and use the sftp command, being sure to supply your key file to access the server:

```
sftp -i path/to/pem/file/keyfile.pem ubuntu@34.209.238.173:path/to/file/on/server/Col0_subsample_SE.fastp.html path/to/local/folder
```
*NOTE* that the IP address will be different (use the same one you used to connect to the server earlier)

Adapter removal is sometimes performed by the sequencing facility (often platform-dependent) or by researchers prior to submission to SRA, so often no adapter sequence may be detected.

### Core Tasks:
* Inspect the `fastp` documentation
* run `fastp` on a single-end read file or paired-end read files
* download the output HTML report to your local machine to view it

## 3. Indexing the genome and aligning reads with Bismark

While multiple alignment tools for BSseq data have been developed (e.g. Biscuit, BWA-meth), Bismark is the best established and most widely used. Briefly, Bismark performs in-silico bisulfite conversion (C->T and G->A) of the input reads and reference genome to facilitate alignment using one of three possible aligners (Bowtie2, Hisat2, or minimap2). High confidence alignments are then compared to the normal genomic sequence (i.e. with Cs included again) to infer the methylation state. By default, Bismark uses the Bowtie2 aligner.

We first need to make an index for the Arabidopsis reference genome. Indexing is a strategy that alignment tools use to speed up alignment and is akin to an index at the end of a book, i.e. allowing specific short sequences to be found rapidly. Read more about genome indexing here: https://pmc.ncbi.nlm.nih.gov/articles/PMC2836519/
Here, we ask Bismark to generate and index a bisulfite-converted genome for use with Bowtie2:
```
# make a directory for the genome: this is the directory we will supply to Bismark
mkdir bisulfite_genome_TAIR10
# copy the genome assembly into the new directory
cp ~/Share/1_Monday/TAIR10_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa bisulfite_genome_TAIR10
# run bismark_genome_preparation
bismark_genome_preparation bisulfite_genome_TAIR10
```
The `bismark_genome_preparation` command tells bismark to make an index of the genome within the 'bisulfite_genome_TAIR10' directory. When it is done, we are ready to align reads to the genome. We can align our quality-trimmed single-end reads as follows:
```
bismark --prefix Col0_align_TAIR10 -p 2 --bam bisulfite_genome_TAIR10 Col0_subsample.trimmed.fastq.gz
```
Or paired-end reads as follows:
```
bismark --prefix Dog4_align_TAIR10 -p 2 --bam bisulfite_genome_TAIR10 \
-1 Dog4_subsample_1.trimmed.fastq.gz -2 Dog4_subsample_2.trimmed.fastq.gz
```
Note that the `--prefix` argument adds a specific string to the beginning of each output file. In any case, the output files are named after the input read files. Adding a prefix can be useful, for example if you want to align the reads to a different assembly (see advanced task, also tomorrow's conversion efficiency task).

When it finishes, you'll see some output files, notably `Col0_align_TAIR10.Col0_subsample.trimmed_bismark_bt2.bam` and `Col0_align_TAIR10.Col0_subsample.trimmed_bismark_bt2_SE_report.txt`.
The .bam file contains the alignments from which methylation calls will be extracted later. The report.txt contains information on the alignment process. You can view it as follows:
```
cat Col0_align_TAIR10.Col0_subsample.trimmed_bismark_bt2_SE_report.txt
```
In particular, the mapping efficiency is derived from the number of reads (or read pairs) with a unique best hit alignment divided by the total number of reads, i.e. the % of reads that could be confidently aligned to the genome. What factors could limit this value?

Finally, let's perform deduplication. This removes reads originating from possible PCR duplicates (inferred due to having identical alignment start or end coordinates) which might otherwise artificially skew subsequent estimates of methylation levels. Note that this may not detect all actual duplicates, given that quality trimming may affect some duplicates differently thus resulting in slightly different alignment start / end coordinates.
```
# remove PCR duplicates from BAM file
deduplicate_bismark --bam --outfile Col0_align_TAIR10 Col0_align_TAIR10.Col0_subsample.trimmed_bismark_bt2.bam
# add --paired if data were paired-end
```

### Core Tasks:
* Examine the `Bismark` documentation: [https://felixkrueger.github.io/Bismark/](https://felixkrueger.github.io/Bismark/)
* Index the Arabidopsis reference genome with `bismark_genome_preparation`
* Align single-end or paired-end reads to the genome with `bismark`
* Examine the alignment report

### Advanced Tasks:
* Copy an alternative genome assembly of one of the non-reference accessions and repeat indexing and alignment. Does the choice of genome affect mapping efficiency?
* Bismark's default behaviour is to try to align reads in end-to-end mode. Try to repeat the alignment step using `--local` mode. How does thia affect mapping efficiency?

Further reading:
https://pmc.ncbi.nlm.nih.gov/articles/PMC4009243/
