#!/bin/bash

# see fastp documentation:
# https://github.com/OpenGene/fastp

N_CORES=2
READ1=Dog4_subsample_1.fastq.gz
READ2=Dog4_subsample_2.fastq.gz
PREFIX=Dog4_subsample

# trim reads with fastp (pair-end)
fastp --thread $N_CORES --in1 $READ1 --out1 ${PREFIX}_1.trimmed.fastq.gz \
    --in2 $READ2 --out2 ${PREFIX}_2.trimmed.fastq.gz \
    -l 35 -h ${PREFIX}_PE.fastp.html \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 10 \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 10 \
    --detect_adapter_for_pe
    
# l = minimum read length to keep
# --cut_front = cut bases from 5' end if below a certain quality
# --cut_tail = cut bases from 3' end if below a certain quality
# --detect_adapter_for_pe = automatically detect adapter sequence (only need to specify for PE reads)
# see also: --cut_right - essentially more aggressive version of --cut_tail
# see also --trim_front1 and --trim_tail1 / --trim_front2 and --trim_tail2 - trim arbitrary number of bases from 5' and 3' ends of read1 / read2