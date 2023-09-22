#!/bin/bash
ref=../bowtie2_ref/mm10/mm10
data=../data
results=../results
threads=16

#Setting up directories
sam_output=$results/alignment/mm10
frag_len_output=$results/alignment/mm10/frag_len

[ ! -d $frag_len_output ] && mkdir -p $frag_len_output

ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

while IFS= read -r sample_name; do
	echo "Getting fragment length..."
	samtools view -F 0x04 $sam_output/${sample_name}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$frag_len_output/${sample_name}_fragmentLen.txt
	echo "Done for $sample_name !"
done < $data/samplelist.txt
