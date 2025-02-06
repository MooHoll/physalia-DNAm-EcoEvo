configfile: "config.yaml"

rule all:
        input:
            'genome/GCF_950023065.1_ihPlaCitr1.1_genomic.fa.sa',
            expand('{sample}_sorted.bam.bai', sample=config["samples"]),
            expand('{sample}.dedup.bam.bai', sample=config["samples"]),
            expand('{sample}_filtered.recode.vcf', sample=config["samples"])
            

rule genome_prep:
    input:
        "genome/GCF_950023065.1_ihPlaCitr1.1_genomic.fa"
    output:
        "genome/GCF_950023065.1_ihPlaCitr1.1_genomic.fa.sa"
    shell:
        """
        bwa index {input}
        """

rule alignment:
    input:
        read1="{sample}_1.fq.gz",
        read2="{sample}_2.fq.gz",
        genome_dir="genome/GCF_950023065.1_ihPlaCitr1.1_genomic.fa",
        genome="genome/GCF_950023065.1_ihPlaCitr1.1_genomic.fa.sa" # needed so the above has run
    output:
        "{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    shell:
        """
        bwa mem -R '{params.rg}' -t 20 {input.genome_dir} {input.read1} {input.read2} \
        | samtools view -Sb - > {output}
        """

# module load bamtools/2.5.2-yu65qvv
# bamtools stats -in *bam

rule sort_bams:
    input:
        "{sample}.bam"
    output:
        "{sample}_sorted.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule index_bams:
    input:
        "{sample}_sorted.bam"
    output:
        "{sample}_sorted.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule rm_duplicates:
    input:
        bam="{sample}_sorted.bam",
        index="{sample}_sorted.bam.bai"
    output:
        bam="{sample}.dedup.bam",
        metrics="{sample}.metrics.txt"
    shell:
        """
        picard \
        MarkDuplicates \
        I={input.bam} \
        O={output.bam} \
        M={output.metrics} \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT
        """

rule index_dedup_bams:
    input:
        "{sample}.dedup.bam"
    output:
        "{sample}.dedup.bam.bai"
    shell:
        """
        samtools index {input}
        """

rule snp_calling_all:
    input:
        bam="{sample}.dedup.bam",
        genome="genome/GCF_950023065.1_ihPlaCitr1.1_genomic.fa"
    output:
        "{sample}.vcf.gz"
    shell:
        """
        freebayes \
        -f {input.genome} \
        -C 2 \
        --min-coverage 5 \
        -u \
        -i \
        -X \
        -b {input.bam} \
        > {output}
        """

rule snp_filtering_all:
    input:
        "{sample}.vcf.gz"
    output:
        "{sample}_filtered.recode.vcf"
    params:
        "{sample}_filtered"
    shell:
        """
        vcftools --gzvcf {input} --max-alleles 2 \
        --minQ 20 --min-meanDP 10 --recode \
        --recode-INFO-all --out {params}
        """