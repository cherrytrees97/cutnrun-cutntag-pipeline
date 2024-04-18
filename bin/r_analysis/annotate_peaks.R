library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(tidyverse)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peak <- readPeakFile("data/pax6_optimal.bed")


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

genomic_feature_plot <- plotAnnoPie(peakAnno)

output_peak_annotation <- as.data.frame(peakAnno)

write.csv(
  output_peak_annotation,
  "results/pax6_optimal_annotated.csv"
)


# Get just the promoters
promoter_1kb <- subset(output_peak_annotation, annotation == "Promoter (<=1kb)")
promoter_1kb_2kb <- subset(output_peak_annotation, annotation == "Promoter (1-2kb)")
promoter_2kb_3kb <- subset(output_peak_annotation, annotation == "Promoter (2-3kb)")

promoter_data <- rbind(promoter_1kb, promoter_1kb_2kb) %>% rbind(promoter_2kb_3kb)

write.csv(
  promoter_data,
  "results/pax6_optimal_promoter_annotated.csv"
)

# Generate a bed file
write.table(
  promoter_data[1:3],
  file = "results/macs2_pax6_peaks_in_promoters.bed",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

# Re-annotating the regions with a Pax6 binding motif
pax6_bound_promoters <- read.table(
  file = "results/meme_denovo_motif/macs_pax6_promoter_only/fimo.tsv",
  header = TRUE
)[3:5]
colnames(pax6_bound_promoters) <- c("seqnames", "start", "end")
promoter_peak <- makeGRangesFromDataFrame(
  pax6_bound_promoters,
  start.field = "start",
  end.field = "end"
)
promoter_peak_anno <- annotatePeak(
  promoter_peak,
  tssRegion=c(-3000, 3000),
  TxDb=txdb,
  annoDb="org.Mm.eg.db"
)
write.csv(
  promoter_peak_anno,
  "results/meme_denovo_motif/macs_pax6_promoter_only/fimo_pax6_bound_peaks_annotated.csv"
)

