# Collect the fragment size information 
library(tidyverse)
library(ggplot2)
library("ggpubr")

sample_file_list = as.list(
  list.files(
    "data/fragment_length",
    pattern = "*.txt",
    full.names = TRUE
  )
)

fragment_lengths = c()

for (file in sample_file_list){
  #Process file names
  file_name = strsplit(file, "/")[[1]][3]
  sample_name = strsplit(file_name, "_")[[1]][1]
  rep_num = strsplit(file_name, "_")[[1]][2]
  
  #Get and format fragment length information
  fragment_lengths = read.table(file, header=FALSE) %>%
    mutate(
      id = paste(sample_name, rep_num, sep = "_"),
      fragLen = V1 %>% as.numeric, 
      fragCount = V2 %>% as.numeric, 
      weight = as.numeric(V2)/sum(as.numeric(V2)), 
      sample = sample_name, 
      replicate = rep_num) %>%
    rbind(fragment_lengths, .)
}

frag_length_violin_plot <- fragment_lengths %>% ggplot(
  mapping = aes(x = id, y = fragLen, weight = weight, fill = sample)) +
  geom_violin(bw=5) + 
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  theme_bw(base_size = 20) +
  ylab("Fragment Length") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

frag_length_graph_plot <- fragment_lengths %>% ggplot(
  mapping = aes(x = fragLen, y = fragCount, group = sample, colour = sample)
) +
  geom_line(size = 1) +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

frag_len_qc_plot <- ggarrange(
  frag_length_violin_plot,
  frag_length_graph_plot,
  ncol = 2,
  nrow = 1,
  common.legend = TRUE, 
  legend="bottom"
)
ggsave(
  "QC_frag_length.png",
  plot = frag_len_qc_plot,
  device = "png",
  path = getwd(),
  units = "in",
  width = 16,
  height = 12,
  bg = "white",
)





