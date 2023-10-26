library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

files <- getSampleFiles()

peak <- readPeakFile("data/Pax6_merged_peaks_bedops.bed")


promoter <- getPromoters(
  TxDb=txdb, 
  upstream=3000, 
  downstream=3000,
  by = "gene"
)

peakAnno <- annotatePeak(
  peak,
  tssRegion=c(-3000, 3000),
  TxDb=txdb,
  annoDb="org.Mm.eg.db"
)

plotAnnoPie(peakAnno)

output_peak_annotation <- as.data.frame(peakAnno)

write.csv(
  output_peak_annotation,
  "Pax6_merged_peaks_beops_annotated.csv"
)