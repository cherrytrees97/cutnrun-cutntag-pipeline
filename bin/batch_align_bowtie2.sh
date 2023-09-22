#!/bin/bash
ref=../bowtie2_ref/mm10/mm10
data=../data
results=../results
threads=16

#Setting up directories
sam_output=$results/alignment/sam
bowtie2_summary=$results/alignment/sam/bowtie2_summary

[ ! -d $sam_output ] && mkdir -p $sam_output
[ ! -d $bowtie2_summary ] && mkdir -p $bowtie2_summary

#Get data
ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

while IFS= read -r sample_name; do
	echo "Aligning $sample_name..."
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${threads} -x ${ref} -1 $data/"$sample_name"_R1_001.fastq.gz -2 $data/"$sample_name"_R2_001.fastq.gz -S $sam_output/"$sample_name"_bowtie2.sam &> $bowtie2_summary/"$sample_name"_bowtie2.txt
	echo "Alignment for $sample_name complete!"
done < $data/samplelist.txt
