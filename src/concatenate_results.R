library(stringr)

# Author: Puja Mehta
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
# Date: November 16, 2023

file_list <- list.files("../results/")
dataset <- read.table(paste0("../results/",file_list[1]), header=TRUE, sep=" ", stringsAsFactors = F)
for (file in file_list[2:length(file_list)]){
  temp_dataset <-read.table(paste0("../results/",file), header=TRUE, sep=" ", stringsAsFactors = F)
  n <- nrow(dataset)+1
  dataset[n,] <- temp_dataset[1,]
}
write.table(dataset,paste0("../results/Results_TwoSampleMR_",Sys.Date(),".tsv"), row.names = F, append = F, quote = F, sep = "\t")
