#!/bin/bash
#batch_frag_len.sh - get counts of fragments for every fragment size for each sample.
#Set these directories first.
data=../data
results=../results

#Setting directory paths and creating directories
input_dir=$results/alignment/
frag_len_output=$results/alignment/frag_len

[ ! -d $frag_len_output ] && mkdir -p $frag_len_output

#Generating sample name list - used for looping.
ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

#Nasty awk function that gives us the fragment length counts
while IFS= read -r sample_name; do
	echo "Getting fragment length..."
	samtools view -F 0x04 $input_dir/${sample_name}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$frag_len_output/${sample_name}_fragmentLen.txt
	echo "Done for $sample_name !"
done < $data/samplelist.txt
