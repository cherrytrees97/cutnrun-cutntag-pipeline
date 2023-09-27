# CUT&TAG/CUT&RUN processing pipeline
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

TODO: write a section on what to look for when reviewing the QC data.

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

This process has been batched using the script `batch_align_bowtie2.sh`.

### 3: Performing QC using R
Two scripts are used to visualize QC metrics about the mapped reads: 
1. `bowtie2_qc_summary.R` - this script summarizes the QC reports given by bowtie2
2. `bowtie2_fragment_length_analysis.R` - this script checks the fragment lengths obtained after CUT&TAG/CUT&RUN sequencing to see if they match what we expect.

Prior to running `bowtie2_fragment_length_analysis.R`, `batch_frag_len.sh` must be run to extract the number of fragments for a given fragment size for each sample. 

### 4: Convert SAM to BAM and BED file format
After performing the QC, the SAM files were converted to both BAM and BED format. As a part of this process, the reads were also filtered to remove unmapped read pairs, and read pairs that are under 1000 bp. 

This process has been batched in `batch_sam_to_bed.sh`. The core commands are as follows: 
```
#Filter and keep mapped read pairs
samtools view -bS -F 0x04 $input_dir/${sample_name}_bowtie2.sam >$bam_output/${sample_name}_bowtie2.mapped.bam

#Convert into bed file format
bedtools bamtobed -i $bam_output/${sample_name}_bowtie2.mapped.bam -bedpe >$bed_output/${sample_name}_bowtie2.bed

#Keep read pairs on same chromosome and fragment length less than 1000 bp
awk '$1==$4 && $6-$2 < 1000 {print $0}' $bed_output/${sample_name}_bowtie2.bed >$bed_output/${sample_name}_bowtie2.clean.bed

#Only extract fragment related columns
cut -f 1,2,6 $bed_output/${sample_name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$bed_output/${sample_name}_bowtie2.fragments.bed
```

TODO: re-evaluate the filters to see if I can reduce the incidence of false positive peaks...

### 5: Convert to bedgraph format for SEACR input.
Because SEACR only accepts bedgraph output as input, with the data value being the sequencing depth at that range, the BED files need to be converted. `bedtools` has `genomecov`, which can perform this calculation and produce a bedgraph output for `SEACR`.

```
bedtools genomecov -i $bed_output/${sample_name}_bowtie2.fragments.bed -g $data/mm10.chrom.sizes -bg > $bedgraph_output/${sample_name}_fragments.bedgraph
```

This process has been batched in the script `batch_bed_to_bedgraph.sh`. `mm10.chrom.sizes` can be downloaded from [here](https://github.com/igvteam/igv/tree/master/genomes/sizes).