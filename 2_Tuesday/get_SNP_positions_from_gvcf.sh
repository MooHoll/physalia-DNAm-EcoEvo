#!/bin/bash

# get positions of all C-T and G-A SNPs in a sinlge-sample GVCF file produced by GATK Haplotypecaller
# This is a little crude as it only considers biallelic SNP sites, but this simple grep / awk method should retrieve the majority of otherwise confounding SNP sites
zcat Col-0.g.vcf.gz | grep -v "#" | awk '$4=="C"&&$5=="T,<NON_REF>"' | cut -f1,2 > Col-0_SNPs_CtoT.txt
zcat Col-0.g.vcf.gz | grep -v "#" | awk '$4=="G"&&$5=="A,<NON_REF>"' | cut -f1,2 > Col-0_SNPs_GtoA.txt