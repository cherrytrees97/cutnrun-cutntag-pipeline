#!/bin/bash
#4_bed_to_bedgraph.sh

data=$1
results=$2
filter=$3

bed_output=$results/alignment/bed
bedgraph_output=$results/bedgraph
[ ! -d $bedgraph_output ] && mkdir -p $bedgraph_output

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt

# -------------------------------------------------------------------------------------------------
# 3: Conversion of BED files to BEDGRAPH files for input into software
# -------------------------------------------------------------------------------------------------
while IFS= read -r sample_name; do
	bedtools genomecov -i $bed_output/${sample_name}_bowtie2.fragments.$filter.bed -g $data/mm10.chrom.sizes -bg > $bedgraph_output/${sample_name}_fragments.bedgraph
done < $data/all_samplelist.txt