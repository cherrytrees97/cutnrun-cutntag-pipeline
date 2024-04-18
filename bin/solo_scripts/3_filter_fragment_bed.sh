#!/bin/bash
#batch_sam_to_bed.sh - convert SAM files to BAM and BED files, with filtering.
#Set these directories first. 
data=$1
results=$2

#Setting directory paths and creating directories
input_dir=$results/alignment
bam_output=$results/alignment/bam
bed_output=$results/alignment/bed

[ ! -d $bam_output ] && mkdir -p $bam_output
[ ! -d $bed_output ] && mkdir -p $bed_output

# -------------------------------------------------------------------------------------------------
# 2: Convert SAM files to BAM and BED files, with filtering.
# -------------------------------------------------------------------------------------------------

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt

while IFS= read -r sample_name; do
	echo "Working on $sample_name"
	#Filter and keep mapped read pairs
    # -F 0x04 removes any unmapped read pairs
	echo "Filtering and converting mapped reads to BAM file format..."
	samtools view -bS -F 0x04 $sam_output/${sample_name}_bowtie2.sam >$bam_output/${sample_name}_bowtie2.mapped.bam
    #Filter blacklisted regions out
    echo "Filtering out blacklisted regions from BAM files.."
    bedtools intersect -v -abam $bam_output/${sample_name}_bowtie2.mapped.bam -b $data/mm10-blacklist.v2.bed > $bam_output/${sample_name}_bowtie2.mapped.blfilter.bam
	#Convert into bed file format
	echo "Converting BAM to BED..."
	bedtools bamtobed -i $bam_output/${sample_name}_bowtie2.mapped.blfilter.bam -bedpe >$bed_output/${sample_name}_bowtie2.bed
	#Keep read pairs on same chromosome and fragment length less than 120 bp
	#120 bp is the original parameter specified in the CUT&RUN paper for TFs
	#The 1000 bp limit previously used was the default for the tutorial using histone mod. data
	echo "Filtering for read pairs on same chromosome and fragment length of less than 121 bp"
	awk '$1==$4 && $6-$2 < 121 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.120.bed
	#Only extract fragment related columns
	echo "Extracting fragment information..."
	cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.120.bed
	#Repeat for the 1000 bp cutof
	echo "Filtering for read pairs on same chromosome and fragment length of less than 121 bp"
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.1000.bed
	#Only extract fragment related columns
	echo "Extracting fragment information..."
	cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.1000.bed

done < $data/all_samplelist.txt

#Create binned BED files for correlation analysis later.
while IFS= read -r sample_name; do
    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bedgraph_output/${sample_name}_bowtie2.fragments.120.bed |
    sort -k1,1V -k2,2n |
    uniq -c |
    awk -v OFS="\t" '{print $2, $3, $1}' | 
    sort -k1,1V -k2,2n  > $bedgraph_output/${sample_name}_bowtie2.fragments120Count.bin$binLen.bed
done < $data/all_samplelist.txt
while IFS= read -r sample_name; do
    awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bedgraph_output/${sample_name}_bowtie2.fragments.1000.bed |
    sort -k1,1V -k2,2n |
    uniq -c |
    awk -v OFS="\t" '{print $2, $3, $1}' | 
    sort -k1,1V -k2,2n  > $bedgraph_output/${sample_name}_bowtie2.fragments1000Count.bin$binLen.bed
done < $data/all_samplelist.txt