library(tidyverse)

setwd('/Users/dwyerdeighan/Desktop/post-processing/chrest')

#all_files=c('CollatedDatasets/perturbed1_collated.csv.gz', 'CollatedDatasets/perturbed2_collated.csv.gz', 'CollatedDatasets/TChem_collated.csv.gz')

all_files = commandArgs(trailingOnly=T) # like: sys.argv
master_df = all_files %>% sort() %>% map2_dfr(., 1:len(.), ~read_csv(.x, lazy=F) |> mutate(time_key=.y)) %>%
    arrange(time, x, y, z) # we do time sorting because it makes more sense to leave last bit as holdout than transience
write.csv(master_df, file=gzfile('ablate_collated.csv.gz'), row.names = F)
#write.csv(master_df, file=gzfile('TChem+Perturbedx2_collated.csv.gz'), row.names = F)
