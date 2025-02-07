configfile: "config.yaml"

rule all:
        input:
            'TAIR10_genome/Bisulfite_Genome',
            'TAIR10_chloroplast/Bisulfite_Genome',
            expand('chloroplast.{sample}_1_bismark_bt2_pe.bam', sample=config["samples"]),
            expand('{sample}.CpG_report.merged_CpG_evidence.cov', sample=config["samples"])

rule genome_prep:
    input:
        "./TAIR10_genome"
    output:
        "TAIR10_genome/Bisulfite_Genome"
    shell:
        """
        bismark_genome_preparation {input}
        """

rule chloroplast_genome_prep:
    input:
        "./TAIR10_chloroplast"
    output:
        "TAIR10_chloroplast/Bisulfite_Genome"
    shell:
        """
        bismark_genome_preparation {input}
        """

rule alignment:
    input:
        read1="{sample}_1.fastq",
        read2="{sample}_2.fastq",
        genome_dir="TAIR10_genome/Bisulfite_Genome" # This is required to make sure previous rules have run
    output:
        "{sample}_1_bismark_bt2_pe.bam"
    shell:
        """
        bismark --multicore 3 --genome TAIR10_genome -1 {input.read1} -2 {input.read2}
        """

rule chloroplast_alignment:
    input:
        read1="{sample}_1.fastq",
        read2="{sample}_2.fastq",
        genome_dir="TAIR10_chloroplast/Bisulfite_Genome" # This is required to make sure previous rules have run
    output:
        "chloroplast.{sample}_1_bismark_bt2_pe.bam"
    shell:
        """
        bismark --multicore 3 --genome TAIR10_chloroplast --prefix chloroplast -1 {input.read1} -2 {input.read2}
        """

rule deduplication:
    input:
        "{sample}_1_bismark_bt2_pe.bam"
    output:
        "{sample}_1_bismark_bt2_pe.deduplicated.bam"
    shell:
        """
        deduplicate_bismark {input}
        """

rule meth_extraction:
    input:
        "{sample}_1_bismark_bt2_pe.deduplicated.bam"
    output:
        "{sample}_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    shell:
        """
        bismark_methylation_extractor \
        --multicore 3 \
        --no_overlap --comprehensive --bedgraph --report --cytosine_report \
        --genome_folder TAIR10_genome/ {input}
        """

rule merge_cpgs:
    input:
        "{sample}_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    output:
        "{sample}.CpG_report.merged_CpG_evidence.cov"
    params:
        "{sample}"
    shell:
        """
        coverage2cytosine \
        -o {params} --merge_CpGs \
        --genome_folder TAIR10_genome/ {input}
        """