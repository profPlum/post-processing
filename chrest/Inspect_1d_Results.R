library(tidyverse)
source('~/.RProfile')

#TC_fn = 'ablate_dummy_batch.csv.gz'
#CT_fn <- 'ablate_control.csv.gz'

#TC_fn = '~/Downloads/Data/5_CPV_data/TChem+CPVs+Zmix_QR_MassR2.csv.gz'
TC_fn = '~/Downloads/Data/10_CPV_data/TChem+CPVs+Zmix_QR_MassR2.csv.gz'
#TC_fn = '~/Downloads/Data/Identity_CPV_data/TChem+CPVs+Zmix_QR.csv.gz'

#CT_fn <- 'CollatedDatasets/ZmixSource=0_NoQR_5CPVs_collated.csv.gz'
CT_fn <- 'CollatedDatasets/ZmixSource=0_10CPVs_collated.csv.gz'
#weights_fn <- '~/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/CTV2_MAE_5CPVs/weights.csv'
sim_type_aliases=list(ctrl='TChem', expr='ChemTab') # for easy comparison across different simulations

setwd('/Users/dwyerdeighan/Desktop/post-processing/chrest')

# don't confuse zmix & mass_CPV_zmix, one is squashed to have normal lm coefficients
CT_df = read_csv_cached(CT_fn) %>% arrange(time) %>% filter(time_key>1)
TC_df = read_csv_cached(TC_fn) %>% select(-starts_with('souspec')) %>%
  arrange(time) %>% filter(time<=max(CT_df$time)&time_key>1)
if ('...1' %in% colnames(TC_df)) TC_df = TC_df %>% select(-...1)
if ('zmix' %in% colnames(TC_df)) {
  lm_df = TC_df %>% select(starts_with('Yi'), zmix)
  Zmix_lm = lm(zmix~.-YiN2-YiAR-YiNO-YiNO2-YiN-YiOH-YiH-YiHO2, data=lm_df)
  print(summary(Zmix_lm))
  stopifnot(summary(Zmix_lm)$r.square.adj>0.999)
  CT_df$zmix=predict(Zmix_lm, CT_df)
}
stopifnot(ncol(CT_df)==ncol(TC_df))

####################### Make Master DF & Post-Process #########################

TC_df$sim_type=sim_type_aliases$ctrl
CT_df$sim_type=sim_type_aliases$expr
master_df = CT_df %>% rbind(TC_df) %>% na.exclude() %>%
  mutate(sim_type=as_factor(sim_type)) %>% 
  rename_all(compose(~sub('_PC_', '_', .x), ~sub('mass_CPV', 'CPV', .x),
                     ~sub('source_CPV', 'dCPV', .x)))
TC_df = master_df %>% filter(sim_type==sim_type_aliases$ctrl)
CT_df = master_df %>% filter(sim_type==sim_type_aliases$expr)
print(colnames(master_df))

# lm(log(abs(souener)+1)~.+.**2, data=cbind(souener=TC_df$souener,lm_df)) %>% summary()

# load('zmix_lm.RData')
# start_TC <- TC_df %>% filter(time_key==2)
# start_CT <- CT_df %>% filter(time_key==2)
# TC_zmix_preds = predict(zmix_lm, newdata=start_TC)
# CT_zmix_preds = predict(zmix_lm, newdata=start_CT)
# stopifnot(R2(TC_df$zmix,TC_zmix_preds)>0.999)
# qqplot(TC_zmix_preds, TC_df$zmix)
# plot(start_TC$zmix, TC_zmix_preds)
# plot(start_TC$zmix)
# plot(TC_zmix_preds)

############################ Plotting Functions ############################

# This also doesn't make sense until we have fixed timesteps (that is beyond the first few time_keys)
start_plot = function(QOI, .time_key=2) master_df %>% filter(time_key==.time_key) %>%
  qplot(x, !!QOI, data=., color=sim_type, main=paste0('start_plot: ', as_label(QOI), ' time_key=', .time_key)) %>%
  print()

qq_on_col = function(col) qqplot(TC_df[[col]], CT_df[[col]], main=paste0('QOI: ', col, ', QQ-plot'),
                                 xlab=paste0(sim_type_aliases$ctrl, '_', col), ylab=paste0(sim_type_aliases$expr,'_', col))

.default_aplha=0.3

QOI_2d_heatmap=function(df, QOI, alpha=.default_aplha) {
  p = ggplot(df, aes(x=x, y=time)) + geom_point(aes(color=!!QOI), alpha=alpha) +
    ggtitle(paste0('QOI: ', as_label(QOI), ', X-pos vs Time, heatmap for 1d Sim')) +
    facet_wrap(~sim_type)
  print(p)
}

# TODO: have times match exactly using either field function or similar
# NOTE: positive residual=overshoot, negative undershoot approx
residual_2d_heatmap = function(df, QOI, relative=F, max_residual=100.0, epsilon=1e-7) {
  CT_QOI = df %>% filter(sim_type==sim_type_aliases$expr) %>% select(!!QOI)
  TC_QOI = df %>% filter(sim_type==sim_type_aliases$ctrl) %>% select(!!QOI)
  
  # Truncate to equal size
  max_rows = pmin(nrow(CT_QOI), nrow(TC_QOI))
  CT_QOI = CT_QOI[1:max_rows,]
  TC_QOI = TC_QOI[1:max_rows,]
  
  prefix=''
  residual = (CT_QOI-TC_QOI) %>% as_vector()
  if (relative) {
    prefix='(Relative) '
    residual = (residual/pmax(as_vector(abs(TC_QOI)), epsilon)) %>% na.exclude() %>% abs()
    residual = residual %>% pmin(max_residual) %>% pmax(-max_residual)
    print(summary(residual))
  }
  df = df %>% filter(sim_type==sim_type_aliases$expr & (1:n()<=length(residual))) %>%
    select(x, time) %>% mutate(residual=residual)
  p = ggplot(df, aes(x=x, y=time)) + geom_point(aes(color=residual), alpha=.default_aplha) +
    ggtitle(paste0(prefix, 'Residual: ', as_label(QOI), ', X-pos vs Time, heatmap for 1d Sim'))
  print(p)
}

bulk_QOI_plots = function(QOIs, qq_plots=T, residual_plots=T, start_plots=T) {
  if (qq_plots) QOIs %>% walk(qq_on_col)
  if (start_plots) QOIs %>% syms %>% walk(~start_plot(.x))
  if (residual_plots) {
    QOIs %>% syms %>% walk(~residual_2d_heatmap(master_df, .x))
    QOIs %>% syms %>% walk(~residual_2d_heatmap(master_df, .x, relative=T))
  }
  
  # heatmap comparisons are awesome!
  QOIs %>% syms %>% walk(~QOI_2d_heatmap(master_df, .x))
}

##########################################################################

QOIs_Yi = c('YiC2H4', 'YiO2', 'YiCO', 'YiCO2', 'YiH2O', 'YiOH', 'YiH2', 'YiCH4')
QOIs=c('zmix', 'souener', 'temp', QOIs_Yi)
(QOIs = QOIs[QOIs %in% colnames(TC_df)])
(QOIs_CPV = master_df %>% select(starts_with('CPV')) %>% summarise_all(sd) %>% 
  as_vector() %>% sort(decreasing=T) %>% head(n=10) %>% names())
(QOIs_CPV_source=sub('CPV_', 'dCPV_', QOIs_CPV))

# # TODO: choose which QOIs & which plots to use
# QOIs %>% walk(qq_on_col)
# bulk_QOI_plots(QOIs)
# bulk_QOI_plots(QOIs_CPV)
# bulk_QOI_plots(QOIs_CPV_source)

################ Investigate Zmix Souener & Temp closely ################

qq_on_col('souener')
qq_on_col('temp')

# look at start plot values to find divergence point
QOIs %>% syms %>% map(start_plot)
start_plot(quo(zmix))

# These are very important plots & summarize the gist of how well the model is doing!
master_df %>% filter(time<1e-5) %>% mutate(souener_log=log(souener-min(souener)+1)) %>%
  QOI_2d_heatmap(quo(souener_log))
master_df %>% filter(time<4e-5) %>% mutate(temp_log=log(temp-min(temp)+1)) %>%
  QOI_2d_heatmap(quo(temp_log))
QOI_2d_heatmap(master_df, quo(temp))
QOI_2d_heatmap(master_df, quo(souener))

master_df %>% filter(time<=2.5e-6) %>% mutate(souener_rank=rank(souener)/n()) %>% QOI_2d_heatmap(quo(souener_rank), alpha=0.6)
master_df %>% filter(time<2.5e-6) %>% mutate(temp_rank=rank(temp)/n()) %>% QOI_2d_heatmap(quo(temp_rank), alpha=0.6)

start_data = master_df %>% filter(time_key==2) %>% arrange(time, x) %>% select(sim_type, starts_with('Yi'))
start_CT = start_data %>% filter(sim_type==sim_type_aliases$expr) %>% select(-sim_type)
start_TC = start_data %>% filter(sim_type==sim_type_aliases$ctrl) %>% select(-sim_type)
all.equal(start_CT, start_TC)

# This won't make sense until we get fixed time steps.
# master_df %>% residual_2d_heatmap(quo(souener), relative=F)
# master_df %>% residual_2d_heatmap(quo(temp), relative=F)

dCPV_TC = TC_df %>% select(starts_with('dCPV'))
dCPV_CT = CT_df %>% select(starts_with('dCPV'))
(xcor_mat = cor(dCPV_TC, dCPV_CT[1:nrow(dCPV_TC),]))
image(xcor_mat, main='Cross correlation between CPV source terms in CT/TC',
      xlab=sim_type_aliases$expr, ylab=sim_type_aliases$ctrl)

############################# Aggregate Plots #############################

# TODO: reuse plot_data to make 3d plot animations! 
agg_QOI_2d_heatmap = function(master_df, QOIs, rank=T, title='QOI Heatmap Comparison') {
  scale_at = switch(rank+1, NULL, QOIs)
  if (rank) title = paste0(title, ' (Ordinal)')
  
  plot_data = master_df %>% #rename_at(QOIs, abbreviate) %>%
    mutate_at(scale_at, base::rank) %>% 
    select(all_of(QOIs), sim_type, time, x, y, z) %>%
    pivot_longer(cols=all_of(QOIs))
  print(summary(plot_data))
  
  # (Inspired) from here: https://stackoverflow.com/questions/68583872/q-q-plot-facet-wrap-with-qq-line-in-r
  p <- ggplot(plot_data, aes(x=x,y=time)) + geom_point(aes(color=value), alpha=.default_aplha) +
    theme(axis.text.x=element_text(angle=45)) + facet_grid(sim_type~name) + 
    ggtitle(title) + labs(color=ifelse(rank, 'value (ordinal)', 'value'))
  ggsave(paste0(gsub(' ', '_', title), '.png'))
  print(p)
}

make_agg_plots = function(master_df) {
  master_df %>% agg_QOI_2d_heatmap(QOIs, rank=T)
  master_df %>% agg_QOI_2d_heatmap(QOIs_Yi, title='Yi Comparison', rank=F)
  master_df %>% agg_QOI_2d_heatmap(QOIs_CPV, title='CPV Comparison', rank=F)
  master_df %>% agg_QOI_2d_heatmap(QOIs_CPV_source, title='CPV_source Comparison', rank=T)
}
master_df %>% filter(time<=1e-5) %>% make_agg_plots()
master_df %>% make_agg_plots()

# parallel coordinate plots, interesting & useful for high-dim data!
n_timesteps=30
CT_df %>% arrange(time) %>% select(starts_with('Yi'), zmix, time) %>%
  group_by(time) %>% summarise_all(mean) %>% head(n=n_timesteps) %>%
  parallel_coordinate_response(quo(time), title=sim_type_aliases$expr, scale=F)
TC_df %>% arrange(time) %>% select(starts_with('Yi'), zmix, time) %>%
   group_by(time) %>% summarise_all(mean) %>% head(n=n_timesteps) %>%
   parallel_coordinate_response(quo(time), title=sim_type_aliases$ctrl, scale=F)

########################## Investigate L1 Constraint ##########################
# 
# weight_matrix = read.csv(weights_fn, row.names=1) %>% as.matrix()
# TC_Yi_matrix = TC_df %>% select(starts_with('Yi')) %>% as.matrix()
# CT_Yi_matrix = CT_df %>% select(starts_with('Yi')) %>% as.matrix()
# CT_CPV_matrix = CT_df %>% select(starts_with('CPV')) %>% as.matrix()
# CT_constrained_CPVs = CT_Yi_matrix%*%weight_matrix
# CT_constrained_CPVs_unity_L2 = (CT_Yi_matrix/apply(CT_Yi_matrix, -2, norm, type='2'))%*%weight_matrix
# 
# mean(apply(CT_Yi_matrix, -2, norm, type='2'))
# 
# CT_constrained_residuals = abs(CT_CPV_matrix-CT_constrained_CPVs)/abs(CT_constrained_CPVs)
# image(CT_constrained_residuals, main='constrained_residuals')
# image(CT_constrained_CPVs, main='constrained_CPVs')
# image(CT_CPV_matrix, main='CPV_matrix')
# View(CT_CPV_matrix)
# View(CT_constrained_CPVs)
# 
# (Yi_l2 = apply(CT_Yi_matrix, -2, norm, type='1'))
# cat('head(Yi_l2):', head(Yi_l2))
# (CPV_l2 = apply(CT_CPV_matrix, -2, norm, type='1'))
# m=lm(Yi_l2~., as_tibble(CT_CPV_matrix))
# summary(m)
# sanity <- abs(CT_constrained_CPVs)<=Yi_l2
# stopifnot(all(sanity))
# stopifnot(all(apply(CT_constrained_CPVs[,-1], -2, norm, type='2')<=Yi_l2))
# # this one is also implied by our proof
# 
# CT_constrained_CPVs_unity_L2 %>% apply(-1, sd)
# CT_constrained_CPVs %>% apply(-1, sd)
# 
# (sanity1 = sweep(head(abs(CT_constrained_CPVs)), 1, head(Yi_l2), FUN='<='))
# (sanity2 = apply(head(abs(CT_constrained_CPVs)), 2, head(Yi_l2), FUN='<=')) # == sanity!
# stopifnot(all(head(sanity)==sanity1) && all(head(sanity)==sanity2))
# 
# # SUPER INTERESTING! LOOK INTO THIS MORE!!
# CT_constrained_CPVs_unity_L2_norm = CT_constrained_CPVs_unity_L2[,-1] %>%
#   apply(-2, norm, type='2')
# hist(CT_constrained_CPVs_unity_L2_norm)
# 
# # # Illustration of Difference between MARGIN in apply & sweep
# # (A = matrix(1:9, 3))
# # (B=sweep(A, MARGIN=1, 1:3, FUN='-'))
# # (C=apply(A, MARGIN=2, 1:3, FUN='-'))
# # stopifnot(all(B==C))
# # # ^ no error they are equal!

# ################# LOESS does exactly the kind of time/space interpolation you want!! #################
# 
# # verified to work 10/12/23
# loess_fit = function(X_df, Y, dims=c('time', 'x', 'y', 'z')) {
#   total_df = as.data.frame(X_df) %>% select(Y=Y, dims) %>%
#     na.exclude() %>% select_if(~sd(.x)>0)
#   print(summary(total_df))
#   loess(Y~., data=total_df)
# }
# 
# #loess_m = TC_df %>% loess_fit('souener')
# 
# fit_field_models = function(df, cols, fit_func=loess_fit) {
#   field_models = cols %>% map(~loess_fit(df,.x))
#   names(field_models)=cols
#   class(field_models)='field_models'
# }
# # field_models = fit_field_models(master_df, cols=QOIs)
# 
# predict.field_models = function(models, newdata) {
#   preds = models %>% map(~predict(.x, newdata=newdata)) %>% reduce(cbind)
#   colnames(preds) = names(models)
# }
# ############################## Error Analysis ##############################

# NOTE: This won't make sense until fixed timesteps are happening

# prefix=''
# SDs = master_df %>% select(starts_with('Yi')) %>% summarise_all(sd) %>% as_vector() %>% sort() %>% .[-(1:5)]
# comp_data = master_df %>% arrange(time, x) %>% select(starts_with('Yi') & any_of(QOIs), time, x) %>% #select_if(~sd(.x)>1e-7 & mean(.x)>0.01) %>%
#   mutate(sim_type=master_df$sim_type)
# mass_frac_data = comp_data %>% filter(sim_type==sim_type_aliases$ctrl) %>% select(-sim_type) %>% arrange(time, x)
# approximation = comp_data %>% filter(sim_type==sim_type_aliases$expr) %>% select(-sim_type) %>% arrange(time, x)
# ALL_MAPE = abs(approximation-(mass_frac_data+1e-8))/(mass_frac_data+1e-8)
# MAPE = apply(ALL_MAPE, -1, mean, na.rm=T)
# MSE = apply((approximation-mass_frac_data)**2, -1, mean)
# MAE = apply(abs(approximation-mass_frac_data), -1, mean)
# VAR = apply(mass_frac_data, -1, var)
# R2 = 1-MSE/VAR
# dotchart(sort(MAE), main=paste(prefix, 'MAE of Surrogate'))
# dotchart(sort(R2), main=paste(prefix, 'R^2 of Surrogate'))
# dotchart(sort(MAPE), main=paste(prefix, 'MAPE of Surrogate'))
# dotchart(sort(VAR**0.5), main=paste(prefix, 'SD of Yis'))
# 
# metric_df = tibble(Species=colnames(mass_frac_data), SD=VAR**0.5, MAE, R2, MAPE)
# plot(metric_df)
# metric_df


