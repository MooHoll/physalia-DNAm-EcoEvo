# DNA Methylation in Ecology and Evolution

### Physalia-Courses 

# Long-reads workshop  // 12th February, 2025

This is the Readme file for the session on Wednesday, where we will cover the ONT analysis. All the data you'll need for this will be stored in: 

`~/Share/3_Wednesday/`

## Overview

We will cover from pod5 modification basecalling to basic visualization of the data, aiming to obtain phased haplotypes in the mealybug Planococcus citri. 

## Dorado base modification basecalling

In case you'd need to download Dorado (it is already installed in the instance for the course), you can rapidly do so by using (this is for Linux, other options on github): 

```
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
tar -xvzf dorado-0.9.1-linux-x64.tar.gz
export LD_LIBRARY_PATH=/pathToDorado/dorado-0.9.1-linux-x64/lib:$LD_LIBRARY_PATH
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

Please *don't run this* as it will just clog the server, but this is how we have done the modification basecalling for the practical:

```
conda activate longreads
~/software/dorado-0.9.1-linux-x64/bin/dorado --reference /home/ubuntu/Share/3_Wednesday/rawdata/GCF_950023065.1_ihPlaCitr1.1_genomic.fa sup,5mCG_5hmCG /home/ubuntu/Share/3_Wednesday/pod5_subset > F4_calls.bam
samtools sort -o F4_sorted.bam F4_calls.bam
samtools index F4_sorted.bam
```

If you have enough CPUs, you can always take advantage of that, using extra cores for samtools with `-@ 10`. 






## Program

### Monday: DNA methylation analysis with short reads (Illumina WGBS and EM-Seq) – 2-7pm Berlin time

1. Presentation: Origin and evolution of DNA methylation in eukaryotes: from parasite defence to gene regulation.
2. Presentation: Introduction to short-read format and processing (WGBS/EM-Seq).
3. Lab 1: Short-read processing:
* QC, filtering, trimming
* Genome preparation
* Alignment, deduplication
