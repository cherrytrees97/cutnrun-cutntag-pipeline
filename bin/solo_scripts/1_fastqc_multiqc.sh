#!/bin/bash
# 1_fastqc_multiqc.sh - generate QC files for review
data=$1
results=$2

fastqc_output=$results/fastqc

[ ! -d $fastqc_output ] && mkdir -p $fastqc_output

fastqc -t 16 $data/*.fastq.gz -o $fastqc_output
multiqc $fastqc_output
mv multiqc* $fastqc_output
