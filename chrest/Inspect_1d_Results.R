library(tidyverse)
library(gganimate)
source('~/.RProfile')

SS_resume_comparison = F

#expr_fn <- 'Perturbed+_FW_collated2/TChem+CPVs+Zmix_MassR2.csv.gz'
#ctrl_fn <- 'Perturbed+_FW_collated1/TChem+CPVs+Zmix_MassR2.csv.gz'

#ctrl_fn = '~/Downloads/Data/5_CPV_data/TChem+CPVs+Zmix_MassR2.csv.gz'
ctrl_fn = '~/Downloads/Data/10_CPV_data/TChem+CPVs+Zmix_MassR2.csv.gz'
expr_fn <- 'CollatedDatasets/ZmixSource=0_NoQR_10CPVs_BUGFIX.csv.gz'

#ctrl_fn = 'ablate_methane_steadystate_collated/TChem+CPVs+Zmix_MassR2.csv.gz'
#expr_fn = 'CollatedDatasets/CT_from_TC_SS_collated.csv.gz'
#SS_resume_comparison = T

weights_fn <- '~/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/CTV2_ZmixSource=0_NoQR_10CPVs/weights.csv'

sim_type_aliases=list(ctrl='TChem', expr='Chemtab') # for easy comparison across different simulations

if (SS_resume_comparison) {
  weights_fn <- 'ablate_methane_steadystate_collated/rot_MassR2.csv.gz'
  sim_type_aliases$expr='Chemtab_Resume'
}
####################### Load & Post-Process Data #########################

setwd('/Users/dwyerdeighan/Desktop/post-processing/chrest')

# don't confuse zmix & mass_CPV_zmix, one is squashed to have normal lm coefficients
Expr_df = read_csv_cached(expr_fn) %>% arrange(time) %>% filter(time_key>1)
Ctrl_df = read_csv_cached(ctrl_fn) %>% select(-starts_with('souspec')) %>%
  arrange(time) %>% filter(SS_resume_comparison | (time<=max(Expr_df$time)&time_key>1))

# Only for SS tests!
if (SS_resume_comparison) {
  Expr_df = Expr_df %>% mutate(time=time+max(Ctrl_df$time)-min(time),
                               time_key=time_key+max(Ctrl_df$time_key)-min(time_key))
}

if ('...1' %in% colnames(Ctrl_df)) Ctrl_df = Ctrl_df %>% select(-...1)
if ('zmix' %in% colnames(Ctrl_df)) {
  lm_df = Ctrl_df %>% select(starts_with('Yi'), zmix)
  Zmix_lm = lm(zmix~.-YiN2-YiAR-YiNO-YiNO2-YiN-YiOH-YiH-YiHO2, data=lm_df)
  print(summary(Zmix_lm))
  stopifnot(summary(Zmix_lm)$r.square.adj>0.999)
  Expr_df$zmix=predict(Zmix_lm, Expr_df)
}
valid_cols = intersect(colnames(Expr_df), colnames(Ctrl_df))
Expr_df = Expr_df %>% select(all_of(valid_cols))
Ctrl_df = Ctrl_df %>% select(all_of(valid_cols))
stopifnot(ncol(Expr_df)==ncol(Ctrl_df))

Ctrl_df$sim_type=sim_type_aliases$ctrl
Expr_df$sim_type=sim_type_aliases$expr
master_df = Expr_df %>% rbind(Ctrl_df) %>% na.exclude() %>%
  mutate(sim_type=as_factor(sim_type)) %>% 
  rename_all(compose(~sub('_PC_', '_', .x), ~sub('mass_CPV', 'CPV', .x),
                     ~sub('source_CPV', 'dCPV', .x))) %>% arrange(time,x,y,z)
Ctrl_df = master_df %>% filter(sim_type==sim_type_aliases$ctrl)
Expr_df = master_df %>% filter(sim_type==sim_type_aliases$expr)
print(colnames(master_df))

############################ Plotting Functions ############################

.default_aplha=0.3
# This also doesn't make sense until we have fixed timesteps (that is beyond the first few time_keys)
start_plot = function(QOI, .time_key=2) master_df %>% filter(time_key==.time_key) %>%
  qplot(x, !!QOI, data=., color=sim_type, alpha=.default_aplha,
        main=paste0('start_plot: ', as_label(QOI), ' time_key=', .time_key)) %>%
  print()

qq_on_col = function(col) qqplot(Ctrl_df[[col]], Expr_df[[col]], main=paste0('QOI: ', col, ', QQ-plot'),
                                 xlab=paste0(sim_type_aliases$ctrl, '_', col), ylab=paste0(sim_type_aliases$expr,'_', col))

QOI_2d_heatmap=function(df, QOI, alpha=.default_aplha, sim_type_facet=T) {
  p = ggplot(df, aes(x=x, y=time)) + geom_point(aes(color=!!QOI), alpha=alpha) +
    ggtitle(paste0('QOI: ', as_label(QOI), ', X-pos vs Time, heatmap for 1d Sim'))
  if ('sim_type' %in% colnames(df) && sim_type_facet) p = p + facet_wrap(~sim_type)
  return(p)
}

## TODO: Debug  package installation to view or save animations...
# QOI_1d_anim=function(df, QOI, sim_type_facet=T) {
#   library(gganimate)
#   df$dt = c(diff(df$time), 0)
#   p = ggplot(df, aes(x=x, y=!!QOI)) + geom_line() + enter_fade() + exit_fade() +
#     transition_states(time,  state_length=dt, transition_length=1e-5) +
#     ggtitle(paste0('QOI: ', as_label(QOI), ', X-pos vs QOI, animation for 1d Sim'))
#   if ('sim_type' %in% colnames(df) && sim_type_facet) p = p + facet_wrap(~sim_type)
#   return(p)
# }

# TODO: have times match exactly using either field function or similar
# NOTE: positive residual=overshoot, negative undershoot approx
residual_2d_heatmap = function(df, QOI, relative=F, max_residual=100.0, epsilon=1e-7) {
  Expr_QOI = df %>% filter(sim_type==sim_type_aliases$expr) %>% select(!!QOI)
  Ctrl_QOI = df %>% filter(sim_type==sim_type_aliases$ctrl) %>% select(!!QOI)
  
  # Truncate to equal size
  max_rows = pmin(nrow(Expr_QOI), nrow(Ctrl_QOI))
  Expr_QOI = Expr_QOI[1:max_rows,]
  Ctrl_QOI = Ctrl_QOI[1:max_rows,]
  
  prefix=''
  residual = (Expr_QOI-Ctrl_QOI) %>% as_vector()
  if (relative) {
    prefix='(Relative) '
    residual = (residual/pmax(as_vector(abs(Ctrl_QOI)), epsilon)) %>% na.exclude() %>% abs()
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
  master_df %>% agg_QOI_2d_heatmap(QOIs_CPV, title='CPV Comparison', rank=T)
  master_df %>% agg_QOI_2d_heatmap(QOIs_CPV_source, title='CPV_source Comparison', rank=T)
}

############################### QOIs: ####################################

QOIs_Yi = c('YiC2H4', 'YiO2', 'YiCO', 'YiCO2', 'YiH2O', 'YiOH', 'YiH2', 'YiCH4')
QOIs=c('zmix', 'souener', 'temp', QOIs_Yi)
(QOIs = QOIs[QOIs %in% colnames(Ctrl_df)])
(QOIs_CPV = master_df %>% select(starts_with('CPV')) %>% summarise_all(sd) %>% 
  as_vector() %>% sort(decreasing=T) %>% head(n=10) %>% names())
(QOIs_CPV_source=sub('CPV_', 'dCPV_', QOIs_CPV))

# # TODO: choose which QOIs & which plots to use
# QOIs %>% walk(qq_on_col)
# bulk_QOI_plots(QOIs)
# bulk_QOI_plots(QOIs_CPV)
# bulk_QOI_plots(QOIs_CPV_source)

############################ Inspect SS Resume Data ############################

if (SS_resume_comparison) {
  # Verify Resume ICs are correct
  QOIs %>% syms %>% map(~start_plot(.x, .time_key=max(Ctrl_df$time_key)))
  
  Expr_df %>% agg_QOI_2d_heatmap(QOIs, rank=T)
  
  Expr_time_df = Expr_df %>% filter(x==min(x)) %>% select(time, time_key) %>% .[1:500,]
  dt=diff(Expr_time_df$time)
  time = Expr_time_df$time[-1]
  summary(dt)
  plot(time, dt, type='l', main='ChemTab dt vs time')
  
  # Inspect SS Convergence in Exactly the Same way that tchemMethaneDiffusionFlame.yaml does
  heatmap_matrix = Expr_df %>% select(temp,x,time_key) %>% pivot_wider(values_from=temp, names_from=time_key) %>% arrange(x)
  rownames(heatmap_matrix)=heatmap_matrix$x
  heatmap_matrix = heatmap_matrix %>% select(-x) %>% as.matrix() %>% t()
  View(heatmap_matrix)
  
  heatmap_diff = diff(heatmap_matrix) # finite difference across rows
  temp_convergence_error_l2 = apply(heatmap_diff, -2, norm)
  temp_convergence_error_l2 = tibble(time_key=as.integer(names(temp_convergence_error_l2)), temp_convergence_error_l2)
  temp_convergence_error_l2 = Expr_df %>% select(time, time_key) %>% unique() %>% right_join(temp_convergence_error_l2)
  
  qplot(time, temp_convergence_error_l2, data=temp_convergence_error_l2)
  stop('Done')
}

################ Investigate Zmix Souener & Temp closely ################

qq_on_col('souener')
qq_on_col('temp')

# look at start plot values to find divergence point
rev(QOIs) %>% syms %>% map(~start_plot(.x, .time_key=2))

# These are very important plots & summarize the gist of how well the model is doing!
master_df %>% filter(time<2.5e-6) %>% mutate(souener_log=log(souener-min(souener)+1)) %>%
  QOI_2d_heatmap(quo(souener_log))
master_df %>% filter(time<2.5e-6) %>% mutate(temp_log=log(temp-min(temp)+1)) %>%
  QOI_2d_heatmap(quo(temp_log))
QOI_2d_heatmap(Expr_df, quo(temp))
# QOI_2d_heatmap(Expr_df, quo(souener))
# QOI_2d_heatmap(master_df, quo(dCPV_zmix))
# QOI_2d_heatmap(master_df, quo(CPV_zmix))
# QOI_2d_heatmap(master_df, quo(zmix))

# Ordinal Temp & Souener
master_df %>% filter(time<=2.5e-6) %>% mutate(souener_rank=rank(souener)/n()) %>% QOI_2d_heatmap(quo(souener_rank), alpha=0.6)
master_df %>% filter(time<2.5e-6) %>% mutate(temp_rank=rank(temp)/n()) %>% QOI_2d_heatmap(quo(temp_rank), alpha=0.6)

# Relative error between ICs:
start_data = master_df %>% filter(time_key==2) %>% arrange(time, x) %>% select(sim_type, starts_with('Yi'))
start_Expr = start_data %>% filter(sim_type==sim_type_aliases$expr) %>% select(-sim_type)
start_Ctrl = start_data %>% filter(sim_type==sim_type_aliases$ctrl) %>% select(-sim_type)
all.equal(start_Expr, start_Ctrl)

# This won't make sense until we get fixed time steps.
# master_df %>% residual_2d_heatmap(quo(souener), relative=F)
# master_df %>% residual_2d_heatmap(quo(temp), relative=F)

# dCPV_Ctrl = Ctrl_df %>% select(starts_with('dCPV'))
# dCPV_Expr = Expr_df %>% select(starts_with('dCPV'))
# (xcor_mat = cor(dCPV_Ctrl, dCPV_Expr[1:nrow(dCPV_Ctrl),]))
# title=paste0('Cross correlation between CPV source terms in ',
#              sim_type_aliases$expr,'/',sim_type_aliases$ctrl)
# image(xcor_mat, main=title, xlab=sim_type_aliases$expr, ylab=sim_type_aliases$ctrl)

############################# Make Aggregate Plots #############################

master_df %>% filter(time<=1e-5) %>% make_agg_plots()
master_df %>% make_agg_plots()

# parallel coordinate plots, interesting & useful for high-dim data!
n_timesteps=30
Expr_df %>% arrange(time) %>% select(starts_with('Yi'), zmix, time) %>%
  group_by(time) %>% summarise_all(mean) %>% head(n=n_timesteps) %>%
  parallel_coordinate_response(quo(time), title=sim_type_aliases$expr, scale=F)
Ctrl_df %>% arrange(time) %>% select(starts_with('Yi'), zmix, time) %>%
   group_by(time) %>% summarise_all(mean) %>% head(n=n_timesteps) %>%
   parallel_coordinate_response(quo(time), title=sim_type_aliases$ctrl, scale=F)

# Check if Zmix changes across time: (it should due to transport)
zmix_sd = Expr_df %>% group_by(x) %>% summarise(CPV_zmix=sd(CPV_zmix))
plot(zmix_sd$x, zmix_sd$CPV_zmix)

########################## Investigate L1 Constraint ##########################

weight_matrix = read.csv(weights_fn, row.names=1) %>% as.matrix()
stopifnot(Expr_df$time==sort(Expr_df$time))
stopifnot(Ctrl_df$time==sort(Ctrl_df$time))
#TC_Yi_matrix = Ctrl_df %>% arrange(time,x) %>% select(starts_with('Yi')) %>% as.matrix()
CT_Yi_matrix = Expr_df %>% select(starts_with('Yi')) %>% rename_all(~sub('Yi', '', .x)) %>% as.matrix()
CT_CPV_matrix = Expr_df %>% select(starts_with('CPV')) %>% as.matrix()
CT_constrained_CPVs = CT_Yi_matrix%*%weight_matrix
stopifnot(colnames(CT_Yi_matrix)==rownames(weight_matrix))
CT_constrained_CPVs_unity_L2 = (CT_Yi_matrix/apply(CT_Yi_matrix, -2, norm, type='2'))%*%weight_matrix

CT_constraint_residuals = abs(CT_CPV_matrix-CT_constrained_CPVs) #/abs(CT_constrained_CPVs)
CT_constraint_residuals_across_time = apply(abs(CT_constraint_residuals), -2, mean)
residuals_df = Expr_df %>% select(time, x,y,z) %>% mutate(MAE_residuals=CT_constraint_residuals_across_time)

QOI_2d_heatmap(residuals_df, quo(MAE_residuals))

residuals_df %>% group_by(time) %>% summarise_at(vars(MAE_residuals), mean) %>%
  ggplot(aes(x=time, y=MAE_residuals)) + geom_line() +
  ggtitle('Constraint Residuals vs Time')

# Get dts:
time_key_time = Expr_df %>% group_by(time_key) %>% summarise(time=mean(time))
dts = diff(time_key_time$time)
summary(dts)
dt = mean(dts)

n_steps = 81 # why?
residuals_df %>% filter(time<=n_steps*dt) %>% arrange(desc(time))

(Yi_l2 = apply(CT_Yi_matrix, -2, norm, type='2'))
cat('head(Yi_l2):', head(Yi_l2))
(CPV_l2 = apply(CT_CPV_matrix, -2, norm, type='2'))
sanity <- abs(CT_constrained_CPVs)<=Yi_l2
stopifnot(all(sanity))
stopifnot(all(apply(CT_constrained_CPVs[,-1], -2, norm, type='2')<=Yi_l2))
# this one is also implied by our proof

# SUPER INTERESTING! LOOK INTO THIS MORE!!
CT_constrained_CPVs_unity_L2_norm = CT_constrained_CPVs_unity_L2[,-1] %>%
  apply(-2, norm, type='2')
hist(CT_constrained_CPVs_unity_L2_norm)

##################### Inspect Zmix Lm: #####################

# load('zmix_lm.RData')
# start_TC <- Ctrl_df %>% filter(time_key==2)
# start_CT <- Expr_df %>% filter(time_key==2)
# TC_zmix_preds = predict(zmix_lm, newdata=start_TC)
# CT_zmix_preds = predict(zmix_lm, newdata=start_CT)
# stopifnot(R2(Ctrl_df$zmix,TC_zmix_preds)>0.999)
# qqplot(TC_zmix_preds, Ctrl_df$zmix)
# plot(start_TC$zmix, TC_zmix_preds)
# plot(start_TC$zmix)
# plot(TC_zmix_preds)

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
# #loess_m = Ctrl_df %>% loess_fit('souener')
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


