#!/bin/bash
data=../data
results=../results

bam_output=$results/alignment/index
deeptools_output=$results/deeptools

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt

while IFS= read -r sample_name; do
    bamCoverage -b $bam_output/${sample_name}_bowtie2.mapped.blfilter.sorted.bam \
    -o $deeptools_output/${sample_name}_coverage.bw \
    -bs 1 \
    -p 16 \
    --normalizeUsing CPM \
    --extendReads \
    --maxFragmentLength 120 \
    --verbose
done < $data/all_samplelist.txt