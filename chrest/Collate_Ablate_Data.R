library(tidyverse)

#input_dir = '/Users/dwyerdeighan/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/_1DSampleDiffusionFlame/flowField_mixtureFraction/'
input_dir = commandArgs(trailingOnly=T)[[1]] # like: sys.argv
all_files = input_dir |> list.files('.*\\.csv\\.gz') %>% file.path(input_dir, .)
master_df = all_files |> map_dfr(read_csv)
write.csv(master_df, file='ablate_collated.csv', row.names = F)
