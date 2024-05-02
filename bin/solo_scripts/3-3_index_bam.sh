#!/bin/bash
data=../data
results=../results

bam_output=$results/alignment/bam
index_output=$results/alignment/index

[ ! -d $index_output ] && mkdir -p $index_output

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt

while IFS= read -r sample_name; do
    echo "Sort bam files first..."
    samtools sort $bam_output/${sample_name}_bowtie2.mapped.blfilter.bam \
    -o $index_output/${sample_name}_bowtie2.mapped.blfilter.sorted.bam \
    -@ 16
    echo "Indexing ${sample_name}_bowtie2.mapped.blfilter.bam..."
    samtools index $index_output/${sample_name}_bowtie2.mapped.blfilter.sorted.bam
    echo "Done indexing."
done < $data/all_samplelist.txt