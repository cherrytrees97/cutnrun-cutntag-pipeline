#!/bin/bash
#batch_align_bowtie2.sh - script to facilitate read mapping using bowtie2. 
#Set these directories first.
ref=../bowtie2_ref/
data=../data
results=../results
threads=16

#Setting directory paths and creating directories.
sam_output=$results/alignment/sam
bowtie2_summary=$results/alignment/sam/bowtie2_summary

[ ! -d $sam_output ] && mkdir -p $sam_output
[ ! -d $bowtie2_summary ] && mkdir -p $bowtie2_summary

#Generating sample name list - used for looping.
ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

#Bowtie2 loop
while IFS= read -r sample_name; do
	echo "Aligning $sample_name..."
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${threads} -x ${ref} -1 $data/"$sample_name"_R1_001.fastq.gz -2 $data/"$sample_name"_R2_001.fastq.gz -S $sam_output/"$sample_name"_bowtie2.sam &> $bowtie2_summary/"$sample_name"_bowtie2.txt
	echo "Alignment for $sample_name complete!"
done < $data/samplelist.txt