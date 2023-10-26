sample_file_list = as.list(
  list.files(
    "data/peak_calling",
    pattern = "*.bed",
    full.names = TRUE
  )
)

peakN = c()
peakWidth = c()

for(file in sample_file_list){
  #Get parameters for samples
  file_name = strsplit(file, "/")[[1]][3]
  sample_name = strsplit(file_name, "_")[[1]][1]
  rep_num = strsplit(file_name, "_")[[1]][2]
  
  # Import data from bed file
  peakInfo = read.table(
    file, 
    header = FALSE, 
    fill = TRUE
  )  %>% mutate(width = abs(V3-V2))
  
  # Extract desired data
  peakN = data.frame(
    peakN = nrow(peakInfo), 
    sampleName = sample_name, 
    replicate = rep_num
  ) %>% rbind(peakN, .)
  peakWidth = data.frame(
    width = peakInfo$width, 
    sampleName = sample_name, 
    replicate = rep_num
  )  %>% rbind(peakWidth, .)
}
peakN %>% select(sampleName, replicate, peakN)