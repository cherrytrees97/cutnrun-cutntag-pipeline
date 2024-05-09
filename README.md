# CUT&TAG/CUT&RUN processing pipeline
Pipeline currently being used to process CUT&RUN/CUT&TAG data. Adapted from Henikoff lab and Harvard Chan Bioinformatics Core resources.
![Summary of pipeline](https://github.com/cherrytrees97/cutnrun-cutntag-pipeline/blob/main/CUT%26TAG%20Analysis%20Pipeline.png)
This document describes how to run the pipeline for analyzing CUT&RUN/CUT&TAG data. 
The pipeline was initially developed mostly from Henikoff lab's tutorial on analyzing [CUT&TAG](https://yezhengstat.github.io/CUTTag_tutorial/#VIII_Differential_analysis) data, but that tutorial focuses specifically on analyzing histone modifications with the technique. An additional tutorial from the Harvard Chan Bioinformatics Core found [here](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/schedule/links-to-lessons.html) was also consulted for considerations regarding transcription factor analysis. Although the Harvard tutorial was initially written for ChIP-Seq analysis, many of the considerations in ChIP-Seq still apply to CUT&RUN and CUT&TAG, and the tutorial has been updated with footnotes that discuss the application of those techniques as well as ATAC-seq.

## Setting up processing environment
This analysis requires the use of a Unix-based operating system.

The following versions of software have been installed in a conda environment: 
1. FastQC >= 0.11.9
2. MultiQC
3. Bowtie2 >= 2.3.4.3 
4. samtools >= 1.10
5. bedtools >= 2.29.1
6. Picard >= 2.18.29
7. deepTools >= 2.0
8. r-essentials
9. ncbi-datasets
10. bedops >= 2.4

R analysis can be conducted in RStudio with the appropriate Renv setup. In the future, I will attempt to generate the R environment necessary in the same conda environment and add the R scripts into the pipeline shell script. 

In addition to above, the peak calling program used in this pipeline is [SEACR](https://github.com/FredHutch/SEACR), which is a purpose-built peak caller for sparse CUT&RUN/CUT&TAG data. For more information on its methodology, refer to its [paper](https://doi.org/10.1186/s13072-019-0287-4). The SEACR repository will need to be cloned prior to analysis.

## Setting up directory structure and filename santization.
The pipeline directory structure is organized in the following format: 

```
bowtie2_ref/
├─ mm10.1.bt2
├─ ...
├─ e.t.c
data/
├─ R1.fastq
├─ R2.fastq
├─ mm10.chrom.sizes
results/
bin/
README.md
```

`data` contains the paired-end read FASTQ files, along with `mm10.chrom.sizes`, retrieved from [here](https://github.com/igvteam/igv/tree/master/genomes/sizes). All output from every step of the pipeline will be stored in the appropriate subdirectory in `results`. The shell scripts to run the pipeline should be stored in `bin`. The shell scripts use relative pathing from the `bin` directory, and therefore must be exectuted from within `bin`. 

Ideally, the FASTQs should be systematically named by condition and replicate ID. If not, the FASTQ filenames should be sanitized. 

```
#NOPE
poopyopoypopoypaoypy_S1_L0001_R1.fastq
#YES
Pax6_rep1_S1_L0001_R1.fastq
```

## Protocol 
### 0: Downloading the appropriate genome files.
#### Download pre-built bowtie2 indices. 
Pre-built bowtie2 indices for a variety of species are available at this [link](https://benlangmead.github.io/aws-indexes/bowtie).
#### Downloading a genome and building its bowtie2 indices.
If you would like to build the indices yourself, install `ncbi-datasets` into your conda environment. More information regarding the installation and operation of NCBI datasets can be found [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/). Using [GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/) as an example, download the genome using the following bash command: 

```
datasets genome taxon mouse --reference --include genome
```
Afterwards, `bowtie2-build` is used to index the genome.
```
bowtie2-build GCF_000001635.27_GRCm39_genomic.fna GRCm39
```

### 1: FastQC and MultiQC
The following Bash command can be run from the pipeline root directory to perform the QC: 

```
mkdir results/fastqc
fastqc -t 16 data/* -o results/fastqc/
multiqc results/fastqc
mv multiqc* results/fastqc
```



### 2: Alignment using Bowtie2
`bowtie2` is the program that we will be using to alignment the reads. The following command is used to start the aligner: 
```
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${threads} -x ${ref} -1 $data/"$sample_name"_R1_001.fastq.gz -2 $data/"$sample_name"_R2_001.fastq.gz -S $sam_output/"$sample_name"_bowtie2.sam &> $bowtie2_summary/"$sample_name"_bowtie2.txt
```
Flag definitions: 
* -I: minimum fragment length
* -X: maximum fragment length
* -p: number of threads to use in this job
* -x: index filename prefix/location (ignore the .X.bt2)
* -1: #1 mate reads
* -2: #2 mate reads
* -S: sam output directory

This process has been batched using the script `2_align_bowtie2.sh`.

### 3: Performing QC using R
Two scripts are used to visualize QC metrics about the mapped reads: 
1. `bowtie2_qc_summary.R` - this script summarizes the QC reports given by bowtie2
2. `bowtie2_fragment_length_analysis.R` - this script checks the fragment lengths obtained after CUT&TAG/CUT&RUN sequencing to see if they match what we expect.

### 4: Convert SAM to BAM and BED file format
After performing the QC, the SAM files were converted to both BAM and BED format. As a part of this process, the reads were also filtered to remove unmapped read pairs, and read pairs that are under 1000 bp. 

This process has been batched in `3_filter_fragment_bed.sh`. The core commands are as follows: 
```
#Filter and keep mapped read pairs
# -F 0x04 removes any unmapped read pairs
samtools view -bS -F 0x04 $sam_output/${sample_name}_bowtie2.sam >$bam_output/${sample_name}_bowtie2.mapped.bam

#Filter blacklisted regions out
bedtools intersect -v -abam $bam_output/${sample_name}_bowtie2.mapped.bam -b $data/mm10-blacklist.v2.bed > $bam_output/${sample_name}_bowtie2.mapped.blfilter.bam

#Convert into bed file format
bedtools bamtobed -i $bam_output/${sample_name}_bowtie2.mapped.blfilter.bam -bedpe >$bed_output/${sample_name}_bowtie2.bed

#Keep read pairs on same chromosome and fragment length less than 1000 bp for TF
awk '$1==$4 && $6-$2 < 1000 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.bed

#Only extract fragment related columns
cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.bed
```


### 4.5 Assess replicate reproducibility
To determine if the biological replicates have high concordance with each other, a correlation analysis should be conducted between all of the sequenced samples. 

```
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $bedgraph_output/${sample_name}_bowtie2.fragments.bed |
sort -k1,1V -k2,2n |
uniq -c |
awk -v OFS="\t" '{print $2, $3, $1}' | 
sort -k1,1V -k2,2n  > $bedgraph_output/${sample_name}_bowtie2.fragmentsCount.bin$binLen.bed
```

To do this, run `3-2_bin_bed.sh`

### 5: Convert to bedgraph format for SEACR input.
Because SEACR only accepts bedgraph output as input, with the data value being the sequencing depth at that range, the BED files need to be converted. `bedtools` has `genomecov`, which can perform this calculation and produce a bedgraph output for `SEACR`.

```
bedtools genomecov -i $bed_output/${sample_name}_bowtie2.fragments.bed -g $data/mm10.chrom.sizes -bg > $bedgraph_output/${sample_name}_fragments.bedgraph
```

This process has been batched in the script `4_bed_to_bedgraph.sh`. `mm10.chrom.sizes` can be downloaded from [here](https://github.com/igvteam/igv/tree/master/genomes/sizes).

### 6: Run SEACR for peak calling. 
To call the peaks with SEACR, the following command is used: 

```
bash $seacr $bedgraph_output/${sample_name}_fragments.bedgraph $control_bedgraph norm stringent $seacr_output/${sample_name}
```

This retrieves all peaks using the `stringent` parameter for SEACR. The command above uses an control IgG CUT&TAG sample to set the peak calling threshold. The control is normalized to the sample prior to peak calling with the `norm` parameter. 

This process has been batched in the script `5_seacr_peak_call.sh`. 

### 7: Finding the intersection and merge of biological replicate peaks.
As a way of narrowing down the candidate list of peaks due to the inherent noise observed even in CUT&TAG for transcription factors, we can take the intersection of the peaks observed in all biological replicates. `bedops intersect` is used over `bedtools intersect` as `bedtools intersect` yields a BED file containing pair-wise intersection between one BED file to the rest of the BED files. `bedops intersect` operates more simply and only returns the exact intersecting intervals across all biological replicates. Note that `bedops intersect` results in the stripping of all peak information from the SEACR peak call files in the resulting output. Additionally, `bedops intersect` does not permit adjusting any parameters for defining intersections; it simply takes the literal intersection.

```
bedops intersect -i $seacr_output/*.stringent.bed
```
