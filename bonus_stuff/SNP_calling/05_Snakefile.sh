configfile: "config.yaml"

rule all:
        input:
            'snps/all_samples_filtered.recode.vcf'

rule snp_calling_all:
    input:
        bam="bam.list",
        genome="genome/NIES_genomic.fa"
    output:
        "snps/all_samples.vcf.gz"
    shell:
        """
        freebayes \
        -f {input.genome} \
        -C 2 \
        --min-coverage 5 \
        -u \
        -i \
        -X \
        --bam-list {input.bam} \
        > {output}
        """

rule snp_filtering_all:
    input:
        "snps/all_samples.vcf.gz"
    output:
        "snps/all_samples_filtered.recode.vcf"
    shell:
        """
        vcftools --gzvcf {input} --max-alleles 2 \
        --minQ 20 --min-meanDP 10 --recode \
        --recode-INFO-all --out all_samples_filtered
        """
        