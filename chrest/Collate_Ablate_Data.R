library(tidyverse)

setwd('/Users/dwyerdeighan/Desktop/post-processing/chrest')

#all_files=c('CollatedDatasets/perturbed1_collated.csv.gz', 'CollatedDatasets/perturbed2_collated.csv.gz', 'CollatedDatasets/TChem_collated.csv.gz')

# post_glob="*.?????.csv.gz"
# data_dir="~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame*/flowField_chemTab"
# all_files=tail(Sys.glob(file.path(data_dir, post_glob)), n=10)

all_files = commandArgs(trailingOnly=T) # like: sys.argv
df = all_files %>% sort() %>% map2_dfr(., 1:len(.), ~read_csv(.x, lazy=F) |> mutate(time_key=.y)) %>%
    arrange(time, x, y, z) # we do time sorting because it makes more sense to leave last bit as holdout than transience

########## New Post-Processing to Derive dt & density ##########

# Verified to work 1/29/24
get_density = function(df) {
  df = df %>% rename_all(~sub('^mass_', '', .x)) # remove mass_CPV prefix
  density_variables = df %>% select(starts_with('density'))
  regular_variables = df[, sub('density_', '', colnames(density_variables))]
  stopifnot(prod(dim(regular_variables))>0)
  stopifnot(sub('Yi|mass', '', colnames(regular_variables))==sub('density_', '', colnames(density_variables)))
  density_estimates = density_variables/regular_variables
  density_estimates = apply(density_estimates, -2, mean, na.rm=T)
  return(density_estimates)
}

density_estimates = get_density(df)
df = df %>% select(-starts_with('density')) %>% mutate(density = density_estimates)

dts = df %>% arrange(time_key) %>% group_by(time_key) %>% summarise(time=mean(time)) %>% mutate(dt=c(diff(time), 0))
dts = dts[-nrow(dts),] # drop final row since dt is poorly defined
df = df %>% inner_join(dts)
################################################################

write.csv(df, file=gzfile('ablate_collated.csv.gz'), row.names = F)
#write.csv(master_df, file=gzfile('TChem+Perturbedx2_collated.csv.gz'), row.names = F)
