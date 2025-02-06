#!/bin/bash

# see fastp documentation:
# https://github.com/OpenGene/fastp

N_CORES=2
READS=Col0_subsample.fastq.gz
PREFIX=Col0_subsample

# trim reads with fastp (single-end)
fastp --thread $N_CORES --in1 $READS --out1 ${PREFIX}.trimmed.fastq.gz \
    -l 35 -h ${PREFIX}_SE.fastp.html \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 10 \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 10
    
# l = minimum read length to keep
# --cut_front = cut bases from 5' end if below a certain quality
# --cut_tail = cut bases from 3' end if below a certain quality
# see also: --cut_right - essentially more aggressive version of --cut_tail
# see also --trim_front1 and --trim_tail1 - trim arbitrary number of bases from 5' and 3' ends of read1