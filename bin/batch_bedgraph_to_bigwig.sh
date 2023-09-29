#!/bin/bash
data=../data
results=../results

bedgraph_output=$results/bedgraph
bigwig_conv_output=$results/bigwig_conv

[ ! -d $bigwig_conv_output ] && mkdir -p $bigwig_conv_output

while IFS= read -r sample_name; do
    bedGraphToBigWig $bedgraph_output/${sample_name}_fragments.bedgraph \
    $data/mm10.chrom.sizes \
    $bigwig_conv_output/${sample_name}_fragments.bw
done < $data/all_samplelist.txt
