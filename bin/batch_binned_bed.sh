#!/bin/bash
data=../data
results=../results

bedgraph_output=$results/bed

binLen=500

while IFS= read -r sample_name; do
    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bedgraph_output/${sample_name}_bowtie2.fragments.bed |
    sort -k1,1V -k2,2n |
    uniq -c |
    awk -v OFS="\t" '{print $2, $3, $1}' | 
    sort -k1,1V -k2,2n  > $bedgraph_output/${sample_name}_bowtie2.fragmentsCount.bin$binLen.bed
done < $data/all_samplelist.txt
