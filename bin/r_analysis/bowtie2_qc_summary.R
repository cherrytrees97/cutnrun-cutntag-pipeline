##=== R command ===## 
library(tidyverse)
library(ggplot2)
library("ggpubr")

## Path to files
sample_file_list = as.list(
  list.files(
    "data/bowtie2_qc",
    pattern = "*.txt",
    full.names = TRUE
  )
)

sample_file_names = as.list(
  list.files(
    "data/bowtie2_qc",
    pattern = "*.txt"
  )
)

# Collect the alignment results from the bowtie2 alignment summary files
alignment_results <- c()

for (file in sample_file_list){
  # Read data from the file into a table
  align_result = read.table(file, header = FALSE, fill = TRUE)
  # Get percentage alignment
  align_rate = substr(
    align_result$V1[6],
    1,
    nchar(as.character(align_result$V1[6]))-1
  )
  # Generate the data frame containing the results and append with rbind
  file_name = strsplit(file, "data/bowtie2_qc/")[[1]][2]
  sample_name = strsplit(file_name, "_")[[1]][1]
  rep_num = strsplit(file_name, "_")[[1]][2]
  
  alignment_results = data.frame(
    id = paste(sample_name, rep_num, sep = "_"),
    target = sample_name,
    rep = rep_num,
    seq_depth = align_result$V1[1] %>% as.character %>% as.numeric,
    mapped_frag_mm10 = align_result$V1[4] %>% as.character %>% as.numeric,
    alignment_rate = align_rate %>% as.numeric
  ) %>% rbind(alignment_results, .)
}

# Plotting some of the data.
# Sequencing depth per million
seq_depth_plot <- alignment_results %>% ggplot(
  mapping = aes(x = target, y = seq_depth/1000000)
) + 
  geom_boxplot() +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth per Million") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")

aligned_frag <- alignment_results %>% ggplot(
  mapping = aes(x = target, y = mapped_frag_mm10/1000000)
) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments per Million") +
  xlab("") + 
  ggtitle("B. Alignable Fragment (mm10)")

percent_aligned <- alignment_results %>% ggplot(
  mapping = aes(x = target, y = alignment_rate)
) + 
  geom_boxplot() +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (mm10)")

bowtie2_qc_plot <- ggarrange(
  seq_depth_plot,
  aligned_frag,
  percent_aligned,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE, 
  legend="bottom"
)

ggsave(
  "bowtie2_qc_plot.png",
  bowtie2_qc_plot,
  device = "png",
  dpi = 300,
  units = "in",
  width = 10,
  height = 8,
  bg = "white",
)

write.csv(alignment_results, "bowtie2_qc_results.csv", row.names = FALSE)
