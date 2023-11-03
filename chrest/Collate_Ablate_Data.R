library(tidyverse)

all_files = commandArgs(trailingOnly=T) # like: sys.argv
master_df = all_files %>% sort() %>% map2_dfr(., 1:len(.), ~read_csv(.x) |> mutate(time_key=.y)) %>%
    arrange(time, x, y, z) # %>% slice_sample(prop=1) # we are shuffling here b/c deterministic splitting will happen later for multiprocessing safeness
write.csv(master_df, file=gzfile('ablate_collated.csv.gz'), row.names = F)
