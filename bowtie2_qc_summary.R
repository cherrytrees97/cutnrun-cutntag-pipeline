##=== R command ===## 
library(tidyverse)

## Path to files
sample_file_list = as.list(
  list.files(
    "data",
    pattern = "*.txt",
    full.names = TRUE
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
  file_name = strsplit(file, "/")[[1]][2]
  sample_name = strsplit(file_name, "_")[[1]][1]
  rep_num = strsplit(file_name, "_")[[1]][2]
  
  alignment_results = data.frame(
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
  ggtitle("B. Algnable Fragment (mm10)")

percent_aligned <- alignment_results %>% ggplot(
  mapping = aes(x = target, y = alignment_rate)
) + 
  geom_boxplot() +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (mm10)")


