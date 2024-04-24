#!/bin/bash
#2_align_bowtie.sh - script to facilitate read mapping using bowtie2. 

#Set these directories first.
data=$1
results=$2
ref=$3
threads=$4

# -------------------------------------------------------------------------------------------------
# 1: Alignment with Bowtie2
# -------------------------------------------------------------------------------------------------
# Directory setup
sam_output=$results/sam
bowtie2_summary=$results/sam/bowtie2_summary

[ ! -d $sam_output ] && mkdir -p $sam_output
[ ! -d $bowtie2_summary ] && mkdir -p $bowtie2_summary

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt

#Bowtie2 loop
# local is used as I did not do any read trimming.
# store the bowtie2 summary results in a file
while IFS= read -r sample_name; do
	echo "Aligning $sample_name..."
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${threads} -x ${ref} -1 $data/"$sample_name"_R1_001.fastq.gz -2 $data/"$sample_name"_R2_001.fastq.gz -S $sam_output/"$sample_name"_bowtie2.sam &> $bowtie2_summary/"$sample_name"_bowtie2.txt
	echo "Alignment for $sample_name complete!"
done < $data/all_samplelist.txt

#Quick QC
frag_len_output=$results/sam/frag_len
[ ! -d $frag_len_output ] && mkdir -p $frag_len_output
#Get fragment length information for each sample.
while IFS= read -r sample_name; do
	echo "Getting fragment length..."
	samtools view -F 0x04 $sam_output/${sample_name}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$frag_len_output/${sample_name}_fragmentLen.txt
	echo "Done for $sample_name !"
done < $data/all_samplelist.txt