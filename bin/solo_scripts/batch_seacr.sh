#!/bin/bash
#batch_seacr.sh - script to run batched peak calling using SEACR.
#Set these diretories first.
data=../data/
results=../results/
#IMPORTANT: SET CORRECT PATH TO SEACR SHELL SCRIPT
seacr=SEACR/SEACR_1.3.sh

#Setting directory paths and creating directories
bedgraph_output=$results/alignment/bedgraph

#Generating sample name list - used for looping.
ls $data | grep fastq | sed s/_R1_001.fastq.gz//g | sed s/_R2_001.fastq.gz//g | sort | uniq > $data/samplelist.txt

#Bedgraph conversion loop.
while IFS= read -r sample_name; do
	bash $seacr $results/bedgraph/${sample_name}_fragments.bedgraph 0.01 non stringent $sample_name
done < $data/samplelist.txt
