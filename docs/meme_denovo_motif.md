# De-novo motif enrichment of our Pax6 data
Analysis start date: 2023-11-03
Analysis end date: ongoing
## Background
One of the next processing steps after peak calling is to identify if there are any enriched motif sequences in our peak data. This will allow us to generate a putative binding profile our transcription factor of interest, Pax6. Since there is a wealth of Pax6 data that is available from other people's work, we can compare the motif that we find here to the motifs that are found in JASPAR to see if there is concordance with our data to theirs. Note that the motifs generated on JASPAR are likely from ChIP-seq data and not CUT&RUN/CUT&TAG, so there is possibility that we may have different values.
## Methods
Four peak call files were used in this analysis: 
SPRI-0613Pax6-JY-072023_S2_L001_qdlm_peaks.narrowPeak
SPRI-0627Pax6a-JY-072023_S4_L001_qdlm_peaks.narrowPeak
SPRI-0627Pax6b-JY-072023_S5_L001_qdlm_peaks.narrowPeak
SPRI-0710Pax6-JY-072023_S7_L001_qdlm_peaks.narrowPeak

After peaking calling using MACS2 (see macs2_comparison.md for the paramaters that were used), ChIP-R was used to perform merge replicate peaks together. The command to run ChIP-R was: 
```
chipr *.narrowPeak -o pax6
```
This produces three output files: `pax6_all.bed`, `pax6_log.txt`, `pax6_optimal.bed`. `pax6_optimal.bed` contains the merged, QC'd peaks from ChIP-R. 

`bedtools getfasta` was used to retrieve the sequences of each peak in `pax6_optimal.bed` on mm10. The mm10 genome was retrieved from UCSC and decompressed using the following command: 
```
 wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
 pigz -d mm10.fa.gz
```
To run `bedtools getfasta`:
```
bedtools getfasta -fi /data2/reference_data/genomes/mm10.fa -fo bedtools_pax6_optimal.fa
```
To run the peak enrichment using MEME: 
```
meme bedtools_pax6_optimal.fa -o meme_output
```