library(tidyverse)
source('~/.RProfile')

#TC_fn <- 'TChem_collated.csv.gz' # default (contains no CPVs)
#TC_fn = '~/Downloads/Data/Identity_CPV_data/TChem+CPVs+Zmix.csv.gz'
#CT_fn <- 'CTV2-R2-10CPV_collated.csv.gz'

TC_fn = '~/Downloads/Data/10_CPV_data/TChem+CPVs+Zmix_MassR2.csv.gz'
CT_fn <- 'CTV2-MAE-10CPV_collated.csv.gz'

setwd('/Users/dwyerdeighan/Desktop/post-processing/chrest')

CT_df = read_csv(CT_fn) %>% arrange(time) %>% filter(time_key>1)
TC_df = read_csv(TC_fn) %>% select(-...1, -starts_with('souspec')) %>%
  arrange(time) %>% filter(time<=max(CT_df$time)&time_key>1)
lm_df = TC_df %>% select(starts_with('Yi'), zmix)
Zmix_lm = lm(zmix~.-YiN2-YiAR-YiNO-YiNO2-YiN-YiOH-YiH, data=lm_df)
summary(Zmix_lm)
CT_df$zmix=predict(Zmix_lm, CT_df)
stopifnot(ncol(CT_df)==ncol(TC_df))

#################### Investigate L1 Constraint ####################

# weights_fn <- '~/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/CTV2_decomp_MAE_10CPVs/weights.csv'
# weight_matrix = read.csv(weights_fn, row.names=1) %>% as.matrix()
# TC_Yi_matrix = TC_df %>% select(starts_with('Yi')) %>% as.matrix()
# CT_Yi_matrix = CT_df %>% select(starts_with('Yi')) %>% as.matrix()
# CT_CPV_matrix = CT_df %>% select(starts_with('CPV')) %>% as.matrix()
# CT_constrained_CPVs = CT_Yi_matrix%*%weight_matrix
# CT_constrained_CPVs2 = (CT_Yi_matrix/apply(CT_Yi_matrix, -2, norm, type='2'))%*%weight_matrix
# 
# CT_constrained_residuals = abs(CT_CPV_matrix-CT_constrained_CPVs)/abs(CT_constrained_CPVs)
# image(CT_constrained_residuals, main='constrained_residuals')
# image(CT_constrained_CPVs, main='constrained_CPVs')
# image(CT_CPV_matrix, main='CPV_matrix')
# View(CT_CPV_matrix)
# View(CT_constrained_CPVs)
# 
# Yi_l2=apply(CT_Yi_matrix**2, -2, sum)**0.5
# CPV_l2 = apply(CT_CPV_matrix, -2, norm, type='2')
# m=lm(Yi_l2~., as_tibble(CT_CPV_matrix))
# summary(m)
# (sanity = sweep(head(abs(CT_constrained_CPVs)), 1, head(Yi_l2), FUN='<='))
# apply(head(abs(CT_constrained_CPVs)), 2, head(Yi_l2), FUN='<=')
# sanity <- abs(CT_constrained_CPVs)<=Yi_l2
# stopifnot(all(sanity))
# stopifnot(all(apply(CT_constrained_CPVs[,-1], -2, norm, type='2')<=Yi_l2))
# # this one is also implied by our proof
# 
# apply(CT_constrained_CPVs2[,-1], -2, norm, type='2')
# 
# # SUPER INTERESTING! LOOK INTO THIS MORE!!
# hist(apply(CT_constrained_CPVs[,-1], -2, norm, type='2')<=Yi_l2)
#
# # Illustration of Difference between MARGIN in apply & sweep
# (A = matrix(1:9, 3))
# (B=sweep(A, MARGIN=1, 1:3, FUN='-'))
# (C=apply(A, MARGIN=2, 1:3, FUN='-'))
# stopifnot(all(B==C)) 
# # ^ no error they are equal!

#####################################################################

TC_df$sim_type='TChem'
CT_df$sim_type='ChemTab'
master_df = CT_df %>% rbind(TC_df) %>% na.exclude() %>%
  mutate(sim_type=as_factor(sim_type)) %>% 
  rename_all(compose(~sub('_PC_', '_', .x), ~sub('mass_CPV', 'CPV', .x),
                     ~sub('source_CPV', 'dCPV', .x)))
TC_df = master_df %>% filter(sim_type=='TChem')
CT_df = master_df %>% filter(sim_type=='ChemTab')
print(colnames(master_df))

#lm(log(abs(souener)+1)~.+.**2, data=cbind(souener=TC_df$souener,lm_df)) %>% summary()

############################ Plotting Functions ############################

start_plot = function(QOI, .time_key=2) master_df %>% filter(time_key==.time_key) %>%
  qplot(x, !!QOI, data=., color=sim_type, main=paste0('start_plot: ', as_label(QOI))) %>%
  print()

qq_on_col = function(col) qqplot(TC_df[[col]], CT_df[[col]], main=paste0('QOI: ', col, ', QQ-plot'),
                                 xlab=paste0('TChem_', col), ylab=paste0('Chemtab_', col))

QOI_2d_heatmap=function(df, QOI) {
  p = ggplot(df, aes(x=x, y=time)) + geom_point(aes(color=!!QOI), alpha=0.1) +
    ggtitle(paste0('QOI: ', as_label(QOI), ', X-pos vs Time, heatmap for 1d Sim')) +
    facet_wrap(~sim_type)
  print(p)
}

QOI_2d_heatmap(master_df, quo(temp))


# TODO: have times match exactly using either field function or similar
# NOTE: positive residual=overshoot, negative undershoot approx
residual_2d_heatmap = function(df, QOI, relative=F, max_residual=100.0,
                               epsilon=1e-7) {
  CT_QOI = df %>% filter(sim_type=='ChemTab') %>% select(!!QOI)
  TC_QOI = df %>% filter(sim_type=='TChem') %>% select(!!QOI)
  
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
  df = df %>% filter(sim_type=='ChemTab' & (1:n()<=length(residual))) %>%
    select(x, time) %>% mutate(residual=residual)
  p = ggplot(df, aes(x=x, y=time)) + geom_point(aes(color=residual), alpha=0.1) +
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

QOIs_Yi = c('zmix', 'YiO2', 'YiCO', 'YiCO2', 'YiH2O', 'YiOH', 'YiH2', 'YiCH4')
QOIs=c('souener', 'temp', QOIs_Yi)
#(QOIs_CPV=paste0('CPV_', 1:9)) # this makes too many assumptions for identity CPV case!
(QOIs_CPV = master_df %>% select(starts_with('CPV')) %>% summarise_all(sd) %>% 
  as_vector() %>% sort(decreasing=T) %>% .[1:9] %>% names())
(QOIs_CPV_source=sub('CPV_', 'dCPV_', QOIs_CPV))

# # TODO: choose which QOIs & which plots to use
# bulk_QOI_plots(QOIs)
# bulk_QOI_plots(QOIs_CPV)
# bulk_QOI_plots(QOIs_CPV_source)

# These are very important plots & summarize the gist of how well the model is doing!
master_df %>% mutate(souener_log=log(souener-min(souener))) %>% QOI_2d_heatmap(quo(souener_log))
master_df %>% mutate(temp_log=log(temp-min(temp))) %>% QOI_2d_heatmap(quo(temp_log))

master_df %>% residual_2d_heatmap(quo(souener), relative=T)

dCPV_TC = TC_df %>% select(starts_with('dCPV'))
dCPV_CT = CT_df %>% select(starts_with('dCPV'))
(xcor_mat = cor(dCPV_TC, dCPV_CT[1:nrow(dCPV_TC),]))
image(xcor_mat, main='Cross correlation between CPV source terms in CT/TC', xlab='CT', ylab='TC')

############################# Aggregate Plots #############################

# TODO: reuse plot_data to make 3d plot animations! 
agg_QOI_2d_heatmap = function(QOIs, rank=T, title='QOI Heatmap Comparison') {
  scale_at = switch(rank+1, NULL, QOIs)
  if (rank) title = paste0(title, ' (Ordinal)')
  
  plot_data = master_df %>% #rename_at(QOIs, abbreviate) %>%
    mutate_at(scale_at, base::rank) %>% 
    select(all_of(QOIs), sim_type, time, x,y,z) %>%
    pivot_longer(cols=all_of(QOIs))
  print(summary(plot_data))
  
  # (Inspired) from here: https://stackoverflow.com/questions/68583872/q-q-plot-facet-wrap-with-qq-line-in-r
  ggplot(plot_data, aes(x=x,y=time)) + geom_point(aes(color=value), alpha=0.1) +
    theme(axis.text.x=element_text(angle=45)) + facet_grid(sim_type~name) + 
    ggtitle(title) + labs(color=ifelse(rank, 'value (ordinal)', 'value'))
}

agg_QOI_2d_heatmap(QOIs, rank=T)
agg_QOI_2d_heatmap(QOIs_Yi, title='Yi Comparison', rank=F)
agg_QOI_2d_heatmap(QOIs_CPV, title='CPV Comparison', rank=F)
agg_QOI_2d_heatmap(QOIs_CPV_source, title='CPV_source Comparison', rank=T)

# parallel coordinate plots, interesting & useful for high-dim data!
CT_df %>% select(starts_with('Yi'), zmix, time) %>% head(n=5000) %>% parallel_coordinate_response(quo(time), title='Chemtab', scale=F)
TC_df %>% select(starts_with('Yi'), zmix, time) %>% head(n=5000) %>% parallel_coordinate_response(quo(time), title='TChem', scale=F)

################# LOESS does exactly the kind of time/space interpolation you want!! #################

# verified to work 10/12/23
loess_fit = function(X_df, Y, dims=c('time', 'x', 'y', 'z')) {
  total_df = as.data.frame(X_df) %>% select(Y=Y, dims) %>%
    na.exclude() %>% select_if(~sd(.x)>0)
  print(summary(total_df))
  loess(Y~., data=total_df)
}

#loess_m = TC_df %>% loess_fit('souener')

fit_field_models = function(df, cols, fit_func=loess_fit) {
  field_models = cols %>% map(~loess_fit(df,.x))
  names(field_models)=cols
  class(field_models)='field_models'
}
# field_models = fit_field_models(master_df, cols=QOIs)

predict.field_models = function(models, newdata) {
  preds = models %>% map(~predict(.x, newdata=newdata)) %>% reduce(cbind)
  colnames(preds) = names(models)
}
############################## Error Analysis ##############################

prefix=''
SDs = master_df %>% select(starts_with('Yi')) %>% summarise_all(sd) %>% as_vector() %>% sort() %>% .[-(1:5)]
comp_data = master_df %>% select(starts_with('Yi')) %>% select(any_of(QOIs)) %>% #select_if(~sd(.x)>1e-7 & mean(.x)>0.01) %>%
  mutate(sim_type=master_df$sim_type)
mass_frac_data = comp_data %>% filter(sim_type=='TChem') %>% select(-sim_type)
approximation = comp_data %>% filter(sim_type=='ChemTab') %>% select(-sim_type)
ALL_MAPE = abs(approximation-(mass_frac_data+1e-8))/(mass_frac_data+1e-8)
MAPE = apply(ALL_MAPE, -1, mean, na.rm=T)
MSE = apply((approximation-mass_frac_data)**2, -1, mean)
MAE = apply(abs(approximation-mass_frac_data), -1, mean)
VAR = apply(mass_frac_data, -1, var)
R2 = 1-MSE/VAR
dotchart(sort(MAE), main=paste(prefix, 'MAE of Surrogate'))
dotchart(sort(R2), main=paste(prefix, 'R^2 of Surrogate'))
dotchart(sort(MAPE), main=paste(prefix, 'MAPE of Surrogate'))
dotchart(sort(VAR**0.5), main=paste(prefix, 'SD of Yis'))

metric_df = tibble(Species=colnames(mass_frac_data), SD=VAR**0.5, MAE, R2, MAPE)
plot(metric_df)
metric_df


