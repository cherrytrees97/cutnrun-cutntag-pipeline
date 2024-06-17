library(tidyverse)
library(ggplot2)
library(corrplot)

reprod = c()
frag_count = NULL

sample_file_list = as.list(
  list.files(
    "data/correlation_analysis",
    pattern = "*.bed",
    full.names = TRUE
  )
)

for (sample_file in sample_file_list) {
  file_name = strsplit(sample_file, "data/correlation_analysis/")[[1]][2]
  sample_name = strsplit(file_name, "_")[[1]][1]
  rep_num = strsplit(file_name, "_")[[1]][2]
  
  if(is.null(frag_count)){
    frag_count = read.table(sample_file, header = FALSE)
    colnames(frag_count) = c("chrom", "bin", paste0(sample_name, "_", rep_num))
  }else{
    frag_count_temp = read.table(sample_file, header=FALSE)
    colnames(frag_count_temp) = c("chrom", "bin", paste0(sample_name, "_", rep_num))
    frag_count = full_join(frag_count, frag_count_temp, by = c("chrom","bin"))
  }
}

M = cor(frag_count %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")

corrplot(
  M, 
  method = "color", 
  outline = T, 
  addgrid.col = "darkgray", 
  order="hclust", 
  addrect = 3, 
  rect.col = "black", 
  rect.lwd = 3,
  cl.pos = "b", 
  tl.col = "indianred4", 
  tl.cex = 1, 
  cl.cex = 1, 
  addCoef.col = "black", 
  number.digits = 2, 
  number.cex = 1, 
  col = colorRampPalette(c("midnightblue","white","darkred"))(100)
)