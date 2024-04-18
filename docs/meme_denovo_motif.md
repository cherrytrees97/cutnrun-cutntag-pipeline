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

## Results
The MEME output produced was very odd. The motifs were extremely large (greater than the 15 to 20 bp that I expected) and also a very small number of sites contributed to each discovered motif. Therefore, I am not sure how to interpret these results and if they're even reliable. MY guess is that because of how low they are, these results are not reflective of Pax6 binding. However, the question remains: is it a quirk with the algorithm used? Or is it due to the quality of our data and the peaks called? The output for MEME can be found in `results/meme_denovo_motif`.

I switched to using the MEME webserver for STREME and FIMO. For STREME, I performed a similar analysis to what was conducted above using MEME. Reviewing the STREME output for the MACS2 ChIP-R QC'd peaks, the motifs that were identified did not match the one deposited in JASPAR. Furthermore, after scanning the 10038 peaks with FIMO, only 567 peaks were identified to contain the JASPAR Pax6 binding motif. 

Unfortunately, this finding may be indicative of poor quality data.

### Examining only peaks in the promoter or enhancer regions.
Another approach to try is to use STREME only for the peaks found in the promoter region or enhancer regions. To do this, I annotated the peaks in `pax6_optimal.bed` in ChIPSeeker and extracted only the peaks that were annotated to be in the promoter regions.
```
promoter_1kb <- subset(output_peak_annotation, annotation == "Promoter (<=1kb)")
promoter_1kb_2kb <- subset(output_peak_annotation, annotation == "Promoter (1-2kb)")
promoter_2kb_3kb <- subset(output_peak_annotation, annotation == "Promoter (2-3kb)")

promoter_data <- rbind(promoter_1kb, promoter_1kb_2kb) %>% rbind(promoter_2kb_3kb)

write.table(
  promoter_data[1:3],
  file = "results/macs2_pax6_peaks_in_promoters.bed",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
```
Afterwards, I followed a similar protocol to above to get the sequences for each peak and ran them through STREME and FIMO.
```
bedtools getfasta -fi /data2/reference_data/genomes/mm10.fa -fo bedtools_pax6_optimal.fa
```

I then took the peaks identified in FIMO and reannotated tem in ChIPSeeker to see what I get. Unfortunately, I got even less peaks (167) than in the previous method. This leads me to believe that there may be a lot of noise in our data that I need to work on removing. MACS2 may be producing way too many false positive peaks.

## Conclusion
Using the MACS2 pipeline still produces a significant amount of noise. The SEACR pipeline may be the approach with integration of all of the libraries into one file prior to peak calling.
