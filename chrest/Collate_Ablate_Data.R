library(tidyverse)

setwd('/Users/dwyerdeighan/Desktop/post-processing/chrest')

#all_files=c('CollatedDatasets/perturbed1_collated.csv.gz', 'CollatedDatasets/perturbed2_collated.csv.gz', 'CollatedDatasets/TChem_collated.csv.gz')

all_files = commandArgs(trailingOnly=T) # like: sys.argv
df = all_files %>% sort() %>% map2_dfr(., 1:len(.), ~read_csv(.x, lazy=F) |> mutate(time_key=.y)) %>%
    arrange(time, x, y, z) # we do time sorting because it makes more sense to leave last bit as holdout than transience

########## New Post-Processing to Derive dt & density ##########
density_Yis = df %>% select(starts_with('density'))
Yis = df %>% select(starts_with('Yi'))
stopifnot(sub('Yi', '', colnames(Yis))==sub('density_', '', colnames(density_Yis)))
density_estimates = density_Yis/Yis
density_estimates = apply(density_estimates, -2, mean, na.rm=T)
df = df %>% select(-starts_with('density')) %>% mutate(density = density_estimates)

dts = df %>% arrange(time_key) %>% group_by(time_key) %>% summarise(time=mean(time)) %>% mutate(dt=c(diff(time), 0))
dts = dts[-nrow(dts),] # drop final row since dt is poorly defined
df = df %>% inner_join(dts)
################################################################

write.csv(df, file=gzfile('ablate_collated.csv.gz'), row.names = F)
#write.csv(master_df, file=gzfile('TChem+Perturbedx2_collated.csv.gz'), row.names = F)
