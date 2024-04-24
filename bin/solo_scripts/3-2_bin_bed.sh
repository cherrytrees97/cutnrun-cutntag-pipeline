#!/bin/bash
#3-2_bin_bed.sh - convert SAM files to BAM and BED files, with filtering.
#Set these directories first. 
#binLen can be about 500
data=$1
results=$2
binLen=$3

#Setting directory paths and creating directories
bed_output=$results/alignment/bed

#Create binned BED files for correlation analysis later.
while IFS= read -r sample_name; do
    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bed_output/${sample_name}_bowtie2.fragments.120.bed |
    sort -k1,1V -k2,2n |
    uniq -c |
    awk -v OFS="\t" '{print $2, $3, $1}' | 
    sort -k1,1V -k2,2n  > $bed_output/${sample_name}_bowtie2.fragments120Count.bin$binLen.bed
done < $data/all_samplelist.txt
while IFS= read -r sample_name; do
    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bed_output/${sample_name}_bowtie2.fragments.1000.bed |
    sort -k1,1V -k2,2n |
    uniq -c |
    awk -v OFS="\t" '{print $2, $3, $1}' | 
    sort -k1,1V -k2,2n  > $bed_output/${sample_name}_bowtie2.fragments1000Count.bin$binLen.bed
done < $data/all_samplelist.txt