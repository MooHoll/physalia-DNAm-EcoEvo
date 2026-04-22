# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

# Long-reads workshop  // 22nd April, 2026
### Instructor Alex de Mendoza (background help Dr. Richard Heery)

This is the Readme file for the session on Wednesday, where we will cover the ONT analysis. All the data you'll need for this will be stored in: 

`~/Share/3_Wednesday/`

## Overview

We will cover from pod5 modification basecalling to basic visualization of the data, aiming to obtain phased haplotypes in the sea anemone *Nematostella vectensis*. 

## Dorado base modification basecalling

In case you'd need to download Dorado (it is already installed in the instance for the course), you can rapidly do so by using (this is for Linux, other options on github): 

```
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.4.0-linux-x64.tar.gz
tar -xvzf dorado-1.4.0-linux-x64.tar.gz
export LD_LIBRARY_PATH=/home/ubuntu/Share/3_Wednesday/software/dorado-1.4.0-linux-x64/lib:$LD_LIBRARY_PATH
```
and then download the base modification models using 

```
dorado download
```

In the workshop we don't have GPUs and it would be computationally very intensive to do basecalling 25 times on the same dataset, but remember that you will need "cuda" ready to interact with the GPUs. This are installed in HPC or you can download using conda: https://anaconda.org/nvidia/cuda. 

Let's have a look at the various modification models in the Dorado github repo: 
[https://github.com/nanoporetech/dorado
](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models)

In a real case scenario, you'll have to look at what is available for the technology you used to sequence your libraries. Not everything is compatible with everything! Dorado and pod5 are very smart and figure out the details of your run for you, so you don't need to worry too much, but still, some combinations are simply not possible. The newest versions will always be better, but try to be consistent across samples (use same model for all your samples!). 

Please **don't run this** as it will just clog the server, but this is how we have done the modification basecalling for the practical:

```
conda activate /opt/miniconda3/envs/longreads
/home/ubuntu/Share/3_Wednesday/software/dorado-1.4.0-linux-x64/bin/dorado basecaller --reference /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna sup,5mCG_5hmCG /home/ubuntu/Share/3_Wednesday/rawdata/pod5 > Sample.bam
samtools sort -o Sample_sorted.bam Sample.bam
samtools index Sample_sorted.bam
samtools view Sample_sorted.bam "NC_064035.1:4600000-4700000" -o Locus_sorted.bam
```

Besides the GPUs, if you have enough CPUs, you can always take advantage of that, using extra cores for samtools with `-@ 10`. 

This bam file contains your reads aligned to the reference genome, with the special flags **MM** and **ML** with the methylation probabilites / values. Let's copy You can have a look at your reads using:

```
ln -s /home/ubuntu/Share/3_Wednesday/rawdata/F4_sorted.bam
samtools index Locus_sorted.bam
samtools view Locus_sorted.bam | head
```

In case you're curious, our reads belong to *Nematostella vectensis* (starlet sea anemone), and for demonstration purposes it is only the region NC_064035.1:4600000-4700000. 


## Base modification calling with ModKit

Modkit is the latest methylation calling algorithm from Nanopore. Before they had modbam2bed, which was probably designed around Guppy outputs. However, both are still compatible.

Modkis has several functions available, one is a quick way to asses global methylation levels in the sampel. 

```
/home/ubuntu/Share/3_Wednesday/software/dist_modkit_v0.6.1_481e3c9/modkit summary -t 2 Locus_sorted.bam
```
This will give an overview of the bam file. 

But what we really want is **pileup**.

As with any softare, good to see options by doing: 
`/home/ubuntu/Share/3_Wednesday/software/dist_modkit_v0.6.1_481e3c9/modkit pileup`

For our test dataset we want to call methylation only on CpGs, this is an insect. Also, dorado used a CG model, so other contexts wouldn't work. 

```
/home/ubuntu/Share/3_Wednesday/software/dist_modkit_v0.6.1_481e3c9/modkit pileup \
--ref /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna \
--threads 2 \
--cpg \
--prefix mCG \
--combine-strands \
--modified-bases 5mC \
/home/ubuntu/Share/3_Wednesday/rawdata/Locus_sorted.bam \
Locus_modkit.bedMethyl
```
This settings only focus on CpGs and combine information from Watson and Crick strand in a single value, all that makes it easier/compact, but depending on analysis, you might want it. Also ignores "h", which is the hydroxmethylation calling. Probably not a great idea for this organism. 

If you had a plant genome or something with non-CG methylation, you could play with **--motif CHH 0** or **--motif CHG 0** (you would need to use a 5mC_5hmC model in dorado). 

To understand the bedMethyl format, check modkit [github repo](https://github.com/nanoporetech/modkit?tab=readme-ov-file#bedmethyl-column-descriptions).

## Visualizing your data in the genome browser

Most important columns, 5 (coverage) and 11 (methylation fraction). Let's make two files to visualize in the genome browser. 
```
cat Locus_modkit.bedMethyl | awk '{print $1,$2,$3,$11}' OFS="\t" > Locus_modkit.mCG.bedGraph
cat Locus_modkit.bedMethyl | awk '{print $1,$2,$3,$5}' OFS="\t" > Locus_modkit.cov.bedGraph
```
These files could be uploaded directly to the genome browser, but to keep it tidy and quick, let's binarise them into bigwig format, using the UCSC tool [bedGraphToBigWig](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig).

```
/home/ubuntu/Share/3_Wednesday/software/bedGraphToBigWig Locus_modkit.mCG.bedGraph /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna.fai Locus_modkit.mCG.bigwig
/home/ubuntu/Share/3_Wednesday/software/bedGraphToBigWig Locus_modkit.cov.bedGraph /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna.fai Locus_modkit.cov.bigwig
```
Download this data into your own computer, together with fasta file and annotation file. 
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/526/225/GCF_932526225.1_jaNemVect1.1/GCF_932526225.1_jaNemVect1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/526/225/GCF_932526225.1_jaNemVect1.1/GCF_932526225.1_jaNemVect1.1_genomic.gtf.gz
gzip -d GCF_932526225.1_jaNemVect1.1_genomic.fna.gz
samtools faidx GCF_932526225.1_jaNemVect1.1_genomic.fna
```
Then download IGV app from: [https://igv.org/doc/desktop/#DownloadPage/].

## Visualizing the per-read methylation data with MethylArtist

[MethylArtist](https://github.com/adamewing/methylartist) is a handy tool to visualize your reads on a certain genomic region.

To get your data ready, you'll need the bam file, and perhaps an indexed gtf file. To index a gtf file do: 
```
gzip -d GCF_932526225.1_jaNemVect1.1_genomic.gtf.gz
sort -k1,1 -k4,4n GCF_932526225.1_jaNemVect1.1_genomic.gtf > sorted.gtf
bgzip sorted.gtf
tabix -p gff sorted.gtf.bgz
```

Methylartist can be installed via conda.

```
conda activate /opt/miniconda3/envs/longreads

methylartist region \
-i NC_064035.1:4626543-4641394 \
-b Locus_sorted.bam \
-r /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna \
-n CG -m m --panelratios 1,2,1,1 \
-o TestRegion.png
# this option below doesn't work for 
-g /home/user1/Share/3_Wednesday/rawdata/sorted.gtf.bgz
```

Download the test region png file to your computer using `scp`.
`
scp -i c1.pem user1@"54.218.18.36:/home/user1/*png" ~/Downloads/TestRegion.png
`

# Phased methylation calling

## Call variants from your bam file
If your long reads have enough covarage and have enough quality (R10 better than R9), you can try to call SNPs on your bam directly with [Clair3](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file) or bcftools.

Here we will do a toy example: 
```
conda activate /opt/miniconda3/envs/clair3_clean

echo -e "NC_064035.1\t4626543\t4641394" > locus.bed

run_clair3.sh --bam_fn=Locus_sorted.bam \
--ref_fn=/home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna \
--threads=2 --platform="ont" \
--model_path=/opt/miniconda3/envs/clair3_clean/bin/models/ont \
--output=clair3_subset --enable_phasing --include_all_ctgs \
--bed_fn=locus.bed

```

The file `locus.bed` is simply to subset the variant calling to small region of the genome. Deleting this option would run it for all contigs/chromosomes, but that takes a long time, so no need for now.

## Phase your reads using WhatsHap

[Whatshap](https://whatshap.readthedocs.io/en/latest/) is a software designed to phase your long reads using a vcf reference. Install via pip works.

First you need to tag your haplogroups found with clair3

```
conda activate /opt/miniconda3/envs/whatshap-env

samtools view -b Locus_sorted.bam -L locus.bed > locus.bam
samtools index locus.bam

whatshap haplotag clair3_subset/phased_merge_output.vcf.gz Locus_sorted.bam \
--reference /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna \
-o Locus_sorted.locus_haplotagged.bam \
--skip-missing-contigs \
--output-threads 2 \
--ignore-read-groups 


samtools index Locus_sorted.locus_haplotagged.bam

samtools view Locus_sorted.locus_haplotagged.bam | grep "HP:i:1" | wc -l
samtools view Locus_sorted.locus_haplotagged.bam | grep "HP:i:2" | wc -l

```
Now your reads for this locus will have a tag depending on the haplogroup they got assigned, the `HP:i:1` will be haplogroup 1, and `HP:i:2` will be haplogroup 2. 

In a real case scenario, you would simply run this with a vcf file for all contigs without subseting your bam. 

## Use methylartist to visualize your locus of interest

Methylartist can take the haplogroup information to give phased methylation information: 

```
conda activate /opt/miniconda3/envs/longreads

methylartist locus \
-i NC_064035.1:4626543-4641394 \
-b Locus_sorted.locus_haplotagged.bam \
-r /home/ubuntu/Share/3_Wednesday/rawdata/GCF_932526225.1_jaNemVect1.1_genomic.fna \
-n CG -m m --panelratios 1,5,1,2,2 \
-o TestLocusHaplo.png --phased \
-g /home/user1/Share/3_Wednesday/rawdata/sorted.gtf.gz

# for a strange region, this didn't work on my end unless I pasted it in a single line code:
methylartist locus -i NC_088680.1:70791301-70816739 -b F4_sorted.locus_haplotagged.bam -r /home/ubuntu/Share/3_Wednesday/rawdata/GCF_950023065.1_ihPlaCitr1.1_genomic.fa -n CG -m m --panelratios 1,5,1,2,2 -o TestLocus3.png --phased -g /home/user1/Share/3_Wednesday/rawdata/GCF_950023065.1_ihPlaCitr1.1_genomic.gtf.bgz 

```
Download the `TestLocusHaplo.png` to see if this locus has any haplotype methylation difference. 


