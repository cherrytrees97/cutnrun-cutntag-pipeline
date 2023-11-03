# Comparison of MACS2 
## Background
The peak calling software most commonly used in ChIP-seq analysis is MACS2. The tutorial from the Harvard Chan bioinformatics core also uses MACS2, and notes that the default parameters can be specified for CUT&RUN data. Therefore, I decided to try using MACS2 as the peak caller for CUT&TAG data (despite the difference, it should be relatively transferrable?) and see if I got a different number of peaks. 

## Protocol
The in-house E18.5 Pax6 CUT&TAG data was used in this analysis. `process_cutntag_cutnrun.sh` was used to process the entire dataset. The `*.fragments.bed` files were the input into MACS2. `batch_macs2_peakcalls.sh` was used to run MACS2 on all of the samples. The parameters are listed here: 

```
macs2 callpeak -t $results/bed/${sample_name}_bowtie2.fragments.bed \
-f BEDPE \
--outdir $peak_output \
-n $sample_name \
-g mm
```

Three output files were generated for each sample: 
| File | Description |
|---|---|
| *.peaks.narrowPeak | Peak information in BED6+4 format. |
| *.peaks.xls | Peak information in Excel spreadsheet format. |
| *.summits.bed | Peak summit information. |

## Results
The output data is stored in a folder named `macs2_peak_calls`. 
The number of peaks identified with MACS2 differs significantly from the number of peaks identified with SEACR (Table 1).

Table 1: Number of called peaks for each sample using SEACR and MACS2. MACS2 peaks were called using the narrowPeaks parameters.
| Sample | New ID | SEACR Peak # | MACS2 Peak # |
|---|---|---|---|
| SPRI-0508Pax6 | Pax6_rep1 | 19369 | 91 |
| SPRI-0613Pax6 | Pax6_rep2 | 10032 | 7773 |
| SPRI-0613Pos | Pos_rep1 | 5781 | 8960 |
| SPRI-0627Pax6a | Pax6dup_rep1 | 12445 | 14344 |
| SPRI-0627Pax6b | Pax6dup_rep2 | 10801 | 34634 |
| SPRI-0627Pos | Pos_rep2 | 4432 | 6991 |
| SPRI-0710Pax6 | Pax6_rep3 | 13481 | 20774 |
| SPRI-0710Pos | Pos_rep3 | 7120 | 10566 |

Note that the positive control is H3K27me3, which are wider peaks and may not be properly captured using MACS2. These will require separate analysis.
An unexpected finding was the fact that Pax6dup_rep1 and Pax6dup_rep2, which are technical replicates, produce such different peak calls when using MACS2. The reason for this is currently unknown. 
Additionally, Pax6_rep1 the sample from whole cerebella has nearly no peaks called by MACS2.

An additional set of optmizations was conducted using this [pre-print](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v1.full) as the basis. For H3K27ac narrow peaks (which I am going to treat as similar to TFs), they used the following parameters: 

```
macs2 callpeak -t $results/bed/${sample_name}_bowtie2.fragments.bed \
-f BEDPE \
--outdir $peak_output \
-n $sample_name \
-g mm \
-q 1e-5 \
--keep-dup all 
--nolambda
--nomodel
```

Therefore, I am going to test the inclusion of each of these parameters in a stepwise fashion. Each of the parameters has been assigned a letter to make naming the output files easier. q = q-value set. d = --keep-dup. l = --nolambda. m = --nomodel

| Sample | New ID | SEACR Peak # | MACS2 Peak # | MACS2 Peak # D | MACS2 Peak # L | MACS2 Peak # M | MACS2 Peak # DLM | MACS2 Peak # QLM | MACS2 Peak # QDLM |
|---|---|---|---|---|---|---|---|---|---|
| SPRI-0508Pax6 | Pax6_rep1 | 19369 | 91 | 211135 | 255 | 91 | 211948 | 142 | 11282 |
| SPRI-0613Pax6 | Pax6_rep2 | 10032 | 7773 | 283161 | 25042 | 7773 | 308884 | 3213 | 49163 |
| SPRI-0710Pax6 | Pax6_rep3 | 13481 | 20774 | 272508 | 25338 | 20774 | 296092 | 11591 | 67406 |
| SPRI-0627Pax6a | Pax6dup_rep1 | 12445 | 14344 | 264696 | 17039 | 14344 | 277238 | 7171 | 58222 |
| SPRI-0627Pax6b | Pax6dup_rep2 | 10801 | 34634 | 207211 | 54059 | 34634 | 224920 | 12956 | 54649 |
| SPRI-0613Pos | Pos_rep1 | 5781 | 8960 | 17211 | 22236 | 8690 | 31951 | 15519 | 9818 |
| SPRI-0627Pos | Pos_rep2 | 4432 | 6991 | 19775 | 20065 | 6991 | 38594 | 9533 | 16204 |
| SPRI-0710Pos | Pos_rep3 | 7120 | 10566 | 24346 | 26099 | 10566 | 48916 | 12096 | 20413 |

Based on the table above, it appears that the QDLM   set of parameters produced peak call files with relatively similar numbers of called peaks between the three biological replicates. The % standard deviation for L, M, QLM conditions were 54, 59, and 51% respectively. The % standard deviation for D, DLM, and QDLM were all 13%. Therefore, we can conclude that keeping duplicate read had the greatest influence on generating consistent peak calls. The average number of called peaks for DLM was 276793 versus 57360 for QDLM. Since QDLM likely removes a large number of false positive peaks by setting the Q-value threshold to be far lower, I decided to move forwards with the QDLM parameter set.

I used ChIP-R to merge the peaks from the four samples, producing a final list of 10000 peaks, with the following output files: `pax6_all.bed`, `pax6_log.txt`, `pax6_optimal.bed`. The command was: 
```
chipr *.narrowPeak -o pax6
```
This program did not take very long to run. 

I ran it through my annotation R script. The distribution of nearest gene element is about the same between MACS2 and SEACR. 

## 2023-10-26, 2023-11-01
My goal is to compare MACS2 and SEACR. Stuff to check: 
1. How many peaks overlap between the MACS2 and SEACR calls?
2. Is the Pax6 binding motif found in both MACS2 and SEACR peak calls?

Starting with #1, I'll be using `bedops` to run these operations. https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html

First set of commands will use a lenient overlap of 1 bp.

```
bedops --intersect *.narrowPeak > Pax6_intersect_macs2.bed
```
Four samples were included in this analysis: 
SPRI-0613Pax6-JY-072023_S2_L001_qdlm_peaks.narrowPeak
SPRI-0627Pax6a-JY-072023_S4_L001_qdlm_peaks.narrowPeak
SPRI-0627Pax6b-JY-072023_S5_L001_qdlm_peaks.narrowPeak
SPRI-0710Pax6-JY-072023_S7_L001_qdlm_peaks.narrowPeak

For the SEACR peaks, I used a similar command: 
```
bedops --intersect *.stringent.bed* > Pax6_intersect_seacr.bed
```
where *.stringent.bed corresponded to the SEACR output of the same 4 samples. 

To check how many SEACR peaks are present in the MACS2 peak file, we use the following command: 
```
bedops -e 1 macs2/Pax6_intersect_macs2.bed seacr/Pax6_intersect_seacr.bed > Pax6_seacr_in_macs2.bed. 
```

The intersection of the four MACS2 peak calls yielded 7967 peaks. The intersection of the four SEACR peak calls yielded 5619 peaks. The result of the comparison was that 5431 of the SEACR peaks were found in the MACS2 peak file. 

## 2023-11-03
The above results are using a rather crude method of identifying overlapping peaks in multiple replicates. Another approach that I need to consider is the use of ChIP-R, which may provide better resolution and a statistical backbone for my peak calling. Therefore, I'll be working on learning about how it works.

## Other code currently not used
```

bedops -e 1 SPRI-0613Pax6-JY-072023_S2_L001_qdlm_peaks.narrowPeak SPRI-0627Pax6a-JY-072023_S4_L001_qdlm_peaks.narrowPeak SPRI-0627Pax6b-JY-072023_S5_L001_qdlm_peaks.narrowPeak SPRI-0710Pax6-JY-072023_S7_L001_qdlm_peaks.narrowPeak > Pax6_macs2_mergedpeaks_bedops.bed

```
