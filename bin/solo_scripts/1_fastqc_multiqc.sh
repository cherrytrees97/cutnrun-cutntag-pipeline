#!/bin/bash
# 1_fastqc_multiqc.sh - generate QC files for review
data=$1
results=$2

fastqc -t 16 $data/* -o $results/fastqc
multiqc $results/fastqc
mv multiqc* $results/fastqc
