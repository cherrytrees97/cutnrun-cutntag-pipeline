#!/bin/bash
#batch_macs2_peakcall.sh - script to batch peak calling w/ MACS2.
#Set these diretories first.
data=../data/
results=../results/

peak_output=$results/macs2_peak_calls
[ ! -d $peak_output ] && mkdir -p $peak_output

#Bedgraph conversion loop.
while IFS= read -r sample_name; do
	macs2 callpeak -t $results/bed/${sample_name}_bowtie2.fragments.bed \
    -f BEDPE \
    --outdir $peak_output \
    -n $sample_name \
    -g mm
done < $data/all_samplelist.txt
