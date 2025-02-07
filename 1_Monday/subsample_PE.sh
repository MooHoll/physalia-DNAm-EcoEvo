#!/bin/bash

reformat.sh \
in=~/Share/1_Monday/raw_data_for_alignment/Dog-4/SRR3301595_1.fastq.gz \
in2=~/Share/1_Monday/raw_data_for_alignment/Dog-4/SRR3301595_2.fastq.gz \
out=Dog4_subsample_1.fastq.gz out2=Dog4_subsample_2.fastq.gz \
samplerate=0.1
