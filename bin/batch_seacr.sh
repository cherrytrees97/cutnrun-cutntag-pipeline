#!/bin/bash
data=../data/
results=../results/alignment/mm10/
seacr=SEACR/SEACR_1.3.sh
threads=16

bedgraph_output=$results/alignment/mm10/bedgraph

ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

while IFS= read -r sample_name; do
	bash $seacr $results/bedgraph/${sample_name}_fragments.bedgraph 0.01 non stringent $sample_name
done < $data/samplelist.txt
