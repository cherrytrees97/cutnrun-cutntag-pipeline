## CUT&TAG/CUT&RUN processing pipeline
This document lists all of the steps that I used to generate the final result files for the Pax6 CUT&TAG dataset. 
The protocol described in this document is adapted from https://yezhengstat.github.io/CUTTag_tutorial/#VIII_Differential_analysis.

### Setting up processing environment
The following versions of software have been installed in a conda environment: 
1. FastQC >= 0.11.9
2. MultiQC
3. Bowtie2 >= 2.3.4.3 
4. samtools >= 1.10
5. bedtools >= 2.29.1
6. Picard >= 2.18.29
7. deepTools >= 2.0

A conda environment file can be found in the root directory of this project.

### Samples
A total of eight samples have been sequenced as a part of this project. The sample preparation was handled by Joanna Yeung and Yuliya.

### Protocol 
#### 0: Downloading the appropriate genome files.
GRCm39 was downloaded using NCBI datasets. The followingcommand was run: 
```
datasets genome taxon mouse --reference --include genome
```

The genome FASTA was stored in `bowtie2_ref/`. Afterwards, `bowtie2-build` was used to index the genome in preparation for alignment. 
```
bowtie2-build GCF_000001635.27_GRCm39_genomic.fna GRCm39
```

Later on, it was noticed that the reference files for CellRanger are for mm10. Therefore, I downloaded the reference bowtie2 files from the following [link](https://benlangmead.github.io/aws-indexes/bowtie).

#### 1: FastQC and MultiQC
Commands were run from the root directory.
```
fastqc -t 16 data/* -o results/fastqc/
multiqc results/fastqc
mv multiqc* results/fastqc
```

#### 2: Alignment using Bowtie2
`bowtie2` is the program that we will be using to alignment the reads. The following command is used to start the aligner: 
```
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/fastq/${histName}_R1.fastq.gz -2 ${projPath}/fastq/${histName}_R2.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
```
Flag definitions: 
* -I: minimum fragment length
* -X: maximum fragment length
* -p: number of threads to use in this job
* -x: index filename prefix/location (ignore the .X.bt2)
* -1: #1 mate reads
* -2: #2 mate reads
* -S: sam output directory

This process has been batched in the script located at `scripts/batch_align_bowtie2.sh`.

Since no spike-in was used, I did not follow Step 3.1.2 in the Henikoff guide.

#### 3: Performing QC using R
Due to the lack of GUI on my server, I will be conducting these steps on R after transferring the `bowtie2_summary` files to my laptop. The scripts and code used in protocol will be saved in the `scripts` directory. Two scripts are used: 
1. `bowtie2_qc_summary.R` - this script summarizes the QC reports given by bowtie2
2. `bowtie2_fragment_length_analysis.R` - this script checks the fragment lengths obtained after CUT&TAG/CUT&RUN sequencing to see if they match what we expect.

Prior to running `bowtie2_fragment_length_analysis.R`, I ran `batch_sam_to_bed.sh` to extract the number of fragments for a given fragment size for each sample. 

#### 4: Convert SAM to BED file format
After performing the QC, the SAM files were converted to BED files using the following script: `batch_sam_to_bed.sh`. The output of this script was stored under `results/alignment/mm10/bed`, with the intermediate BAM files stored in `results/alignment/mm10/bam`.
