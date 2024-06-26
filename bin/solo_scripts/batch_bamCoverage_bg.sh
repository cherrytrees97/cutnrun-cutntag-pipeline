#!/bin/bash
data=../data
results=../results
length=$1

bam_output=$results/alignment/index
deeptools_output=$results/deeptools

[ ! -d $deeptools_output ] && mkdir -p $deeptools_output

#Generating sample name list - used for looping at all steps.
ls $data | grep "fastq" | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/all_samplelist.txt

while IFS= read -r sample_name; do
    bamCoverage -b $bam_output/${sample_name}_bowtie2.mapped.blfilter.sorted.bam \
    -o $deeptools_output/${sample_name}_coverage.bedgraph \
    -bs 1 \
    -p 16 \
    --normalizeUsing CPM \
    --extendReads \
    --maxFragmentLength $length \
    --outFileFormat bedgraph \
    --verbose
done < $data/all_samplelist.txt