#!/bin/bash
#batch_sam_to_bed.sh - convert SAM files to BAM and BED files, with filtering.
#Set these directories first. 
data=../data/
results=../results

#Setting directory paths and creating directories
input_dir=$results/alignment
bam_output=$results/alignment/bam
bed_output=$results/alignment/bed

[ ! -d $bam_output ] && mkdir -p $bam_output
[ ! -d $bed_output ] && mkdir -p $bed_output

#Generating sample name list - used for looping.
ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

#Conversion loop
while IFS= read -r sample_name; do
	echo "Working on $sample_name"
	#Filter and keep mapped read pairs
	echo "Filtering and converting mapped reads to BAM file format..."
	samtools view -bS -F 0x04 $input_dir/${sample_name}_bowtie2.sam >$bam_output/${sample_name}_bowtie2.mapped.bam
	#Convert into bed file format
	echo "Converting BAM to BED..."
	bedtools bamtobed -i $bam_output/${sample_name}_bowtie2.mapped.bam -bedpe >$bed_output/${sample_name}_bowtie2.bed
	#Keep read pairs on same chromosome and fragment length less than 1000 bp
	echo "Filtering for read pairs on same chromosome and fragment length of less than 1000 bp"
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.bed
	#Only extract fragment related columns
	echo "Extracting fragment information..."
	cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.bed
done < $data/samplelist.txt
