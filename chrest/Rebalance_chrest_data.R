library(tidyverse)

#input_fn = '/Users/dwyerdeighan/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/_1DSampleDiffusionFlame/flowField_mixtureFraction/flowField_mixtureFraction.00001.csv.gz'
input_fn = commandArgs(trailingOnly=T)[[1]] # like: sys.argv
cat('input_fn: ', input_fn, '\\n')

# 1900 degrees Kelvin is the approximate temperature of candle flames apparently
is_flame_predicate=function(temp, YiO2) temp>=1900 & YiO2<0.8

# First we use lazy loading and check entire dataset to find
# the distribution of souener (& an ideal rule for classifying air vs flame)
temp_by_group = read_csv(input_fn, col_select=c(group, temp, YiO2)) %>% group_by(group) %>%
    summarize_all(mean) %>% ungroup %>% mutate(is_flame=is_flame_predicate(temp, YiO2))

n_flame = sum(temp_by_group$is_flame) # we want only has much air as there is flame
cat('(before) mean(temp_by_group$is_flame): ', mean(temp_by_group$is_flame), '\n')
temp_by_group = temp_by_group %>% group_by(is_flame) %>% slice_sample(n=n_flame)
cat('(after) mean(temp_by_group$is_flame): ', mean(temp_by_group$is_flame), '\n')
stopifnot(mean(temp_by_group$is_flame)>=0.5)

##### Reduce memory usage: #####
kept_groups=temp_by_group$group
rm(temp_by_group)
gc()
################################

# lazy load so that we can filter
df = read_csv(input_fn, lazy=T) |> filter(group %in% kept_groups) %>%
  mutate(is_flame=is_flame_predicate(temp, YiO2))
print(c('colnames: ', colnames(df)))

# write back to input file
out_fn=sub('\\.csv.*$', '-rebalanced.csv.gz', input_fn)
cat('output file:', out_fn)
write.csv(df, out_fn)
