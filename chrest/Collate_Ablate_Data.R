library(tidyverse)

#input_dir = '/Users/dwyerdeighan/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/_1DSampleDiffusionFlame/flowField_mixtureFraction/'
#input_dir = commandArgs(trailingOnly=T)[[1]] # like: sys.argv
#all_files = input_dir |> list.files('.*[0-9]{5}\\.csv\\.gz') %>% file.path(input_dir, .)
all_files = commandArgs(trailingOnly=T) # like: sys.argv
master_df = all_files %>% sort() %>% map2_dfr(., 1:len(.), ~read_csv(.x) |> mutate(time_key=.y))
write.csv(master_df, file=gzfile('ablate_collated.csv.gz'), row.names = F)