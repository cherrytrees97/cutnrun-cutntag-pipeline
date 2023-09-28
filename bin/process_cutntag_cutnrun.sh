#!/bin/bash
# -------------------------------------------------------------------------------------------------
# process_cutngtag_cutnrun.sh - script to run the entire pipeline.
# 
# Prior to this, basic QC using FASTQC and MultiQC should be conducted.
# Ensure that the conda environment is activated.
# Clone the SEACR repository () to bin or the same directory wwhere this script is being run.
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# 0: Parameters
# -------------------------------------------------------------------------------------------------
data=../data
results=../results
ref=~/Documents/CUTnTAG_CUTnRUN/mm10_bowtie2_indices/mm10
seacr=SEACR/SEACR_1.3.sh
threads=16

#Use to generate two lists
target_ident=Pax6
pos_ident=Pos
target_list=$data/${target_ident}_samplelist.txt
pos_list=$data/${pos_ident}_samplelist.txt

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt
ls $data | grep "fastq" | grep "${sample_ident}" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $target_list
ls $data | grep "fastq" | grep "${pos_ident}" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $pos_list

# -------------------------------------------------------------------------------------------------
# 1: Alignment with Bowtie2
# -------------------------------------------------------------------------------------------------
# Directory setup
sam_output=$results/sam
bowtie2_summary=$results/sam/bowtie2_summary

[ ! -d $sam_output ] && mkdir -p $sam_output
[ ! -d $bowtie2_summary ] && mkdir -p $bowtie2_summary

#Bowtie2 loop
# local is used as I did not do any read trimming.
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

# -------------------------------------------------------------------------------------------------
# 2: Convert SAM files to BAM and BED files, with filtering.
# -------------------------------------------------------------------------------------------------
bam_output=$results/bam
bed_output=$results/bed

[ ! -d $bam_output ] && mkdir -p $bam_output
[ ! -d $bed_output ] && mkdir -p $bed_output

#Run for TFs
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
	awk '$1==$4 && $6-$2 < 121 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.bed
	#Only extract fragment related columns
	echo "Extracting fragment information..."
	cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.bed
done < $target_list

#Run for positive control
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
	#Keep read pairs on same chromosome and fragment length less than 1000 bp
	echo "Filtering for read pairs on same chromosome and fragment length of less than 1000 bp"
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.bed
	#Only extract fragment related columns
	echo "Extracting fragment information..."
	cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.bed
done < $pos_list

# -------------------------------------------------------------------------------------------------
# 3: Conversion of BED files to BEDGRAPH files for input into software
# -------------------------------------------------------------------------------------------------
bedgraph_output=$results/bedgraph
[ ! -d $bedgraph_output ] && mkdir -p $bedgraph_output

while IFS= read -r sample_name; do
	bedtools genomecov -i $bed_output/${sample_name}_bowtie2.fragments.bed -g $data/mm10.chrom.sizes -bg > $bedgraph_output/${sample_name}_fragments.bedgraph
done < $data/all_samplelist.txt

# -------------------------------------------------------------------------------------------------
# 4: Run the peak caller SEACR
# -------------------------------------------------------------------------------------------------
seacr_output=$results/peak_calls
[ ! -d $seacr_output ] && mkdir -p $seacr_output

#Bedgraph conversion loop.
while IFS= read -r sample_name; do
	bash $seacr $bedgraph_output/${sample_name}_fragments.bedgraph 0.01 non stringent $seacr_output/${sample_name}
done < $data/all_samplelist.txt
