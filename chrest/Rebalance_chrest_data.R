library(tidyverse)

input_fn = commandArgs(trailingOnly=T)[[1]] # like: sys.argv
cat('input_fn: ', input_fn, '\\n')

# is_flame_predicate=function(souener) {
#   #souener_abs_min = max(abs(souener))*0.01 # based on kenny's advice
#   souener_abs_min = exp(log(max(abs(souener)))/2) # my version
#   return(abs(souener)>souener_abs_min)
# }

# loads hot_lm wrapped in: is_flame_predicate() func
load('hot_lm.RData')

# we export the predicate directly for later use
is_flame_predicate=function(souener) {
  preds=predict(hot_lm, newdata=tibble(souener), type='response')
  stopifnot(0.0<preds & preds<1.0)
  return(preds>0.5)
}

# # loads hot_lm wrapped in: is_flame_predicate() func
# load('hot_lm-is_flame_predicate.RData')

# First we use lazy loading and check entire dataset to find
# the distribution of souener (& an ideal rule for classifying air vs flame)
souener_by_group = read_csv(input_fn, col_select=c(group, souener), lazy=T) %>% group_by(group) %>%
  summarize_all(mean) %>% ungroup %>% mutate(is_flame=is_flame_predicate(souener))

n_flame = sum(souener_by_group$is_flame) # we want only has much air as there is flame
cat('mean(souener_by_group$is_flame): ', mean(souener_by_group$is_flame), '\n')
#stopifnot(mean(souener_by_group$is_flame)<0.5)
kept_groups = souener_by_group %>% group_by(is_flame) %>% slice_sample(n=n_flame) %>% ungroup
stopifnot(mean(kept_groups$is_flame)==0.5)

##### Reduce memory usage: ##### 
kept_groups=kept_groups$group
rm(souener_by_group)
gc()
################################

# lazy load so that we can filter
df = read_csv(input_fn) |> filter(group %in% kept_groups) 
#%>% group_by(group) %>% summarize_all(mean) %>% ungroup
print(c('colnames', colnames(df)))

# write back to input file
out_fn=sub('\\.csv.*$', '-rebalanced\\.csv.gz', input_fn)
cat('output file:', out_fn)
write.csv(df, out_fn)
