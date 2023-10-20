#!/bin/bash
data=../data
results=../results

bam_output=$results/bam

while IFS= read -r sample_name; do
    echo "Sort bam files first..."
    samtools sort $bam_output/${sample_name}_bowtie2.mapped.blfilter.bam \
    -o $bam_output/${sample_name}_bowtie2.mapped.blfilter.sorted.bam
    echo "Indexing ${sample_name}_bowtie2.mapped.blfilter.bam..."
    samtools index $bam_output/${sample_name}_bowtie2.mapped.blfilter.sorted.bam
    echo "Done indexing."
done < $data/all_samplelist.txt