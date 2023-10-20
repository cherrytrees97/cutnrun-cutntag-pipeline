#!/bin/bash
#batch_bed_to_bedgraph.sh - convert fragment .bed file to bedgraph format for SEACR.
#Set these directories first.
data=../data/
results=../results

#Setting directory paths and creating directories
bed_output=$results/alignment/bed
bedgraph_output=$results/alignment/bedgraph

[ ! -d $bedgraph_output ] && mkdir -p $bedgraph_output

#Generating sample name list - used for looping.
ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

#Bedgraph conversion loop
while IFS= read -r sample_name; do
	bedtools genomecov -i $bed_output/${sample_name}_bowtie2.fragments.bed -g $data/mm10.chrom.sizes -bg > $bedgraph_output/${sample_name}_fragments.bedgraph
done < $data/samplelist.txt
