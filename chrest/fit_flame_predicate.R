library(tidyverse)

df = read_csv(file.choose())
summary(df)
colnames(df)

################# Plot Important Yis #################
Yi_importance = df %>% select(starts_with('Yi')) %>% summarize_all(sd) %>%
  sort() %>% as_vector() %>% tail(n=10)
names(Yi_importance)=names(Yi_importance) %>% sub('Yi', '', .) %>% abbreviate()
dotchart(Yi_importance)
#####################################################

# # Considering using regular linear model?
# T_lm = lm(T~log(abs(souener)+1)+YiO2+souspecO2+Zmix, data=df)
# summary(T_lm)

# hot_lm is most straight-forward using it for now
hot_lm = glm(T>mean(T)~log(abs(souener)+1), data=df)
summary(hot_lm)

# we export the predicate directly for later use
is_flame_predicate=function(souener) {
  preds=predict(hot_lm, newdata=tibble(souener), type='response')
  stopifnot(0.0<preds & preds<1.0)
  return(preds>0.5)
}

save(hot_lm, file='hot_lm.RData')
save(is_flame_predicate, file='hot_lm-is_flame_predicate.RData')
# ^ 2nd is maybe more advanced pattern?
