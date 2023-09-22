#!/bin/bash
ref=../bowtie2_ref/mm10/mm10
data=../data/
results=../results
threads=16

bam_output=$results/alignment/mm10/bam
bed_output=$results/alignment/mm10/bed
bedgraph_output=$results/alignment/mm10/bedgraph

[ ! -d $bam_output ] && mkdir -p $bam_output
[ ! -d $bed_output ] && mkdir -p $bed_output
[ ! -d $bedgraph_output ] && mkdir -p $bedgraph_output

ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

while IFS= read -r sample_name; do
	bedtools genomecov -i $bed_output/${sample_name}_bowtie2.fragments.bed -g $data/mm10.chrom.sizes -bg > $bedgraph_output/${sample_name}_fragments.bedgraph
done < $data/samplelist.txt

mv $bed_output/*.bedgraph $bedraph_output
