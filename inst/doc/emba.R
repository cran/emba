## ----setup, message=FALSE-----------------------------------------------------
# libraries
library(emba)
library(usefun)
library(dplyr)
library(knitr)
library(Ckmeans.1d.dp)

# wrapper to change printing to invisible
pr = function(x) invisible(x)

## ----input-1------------------------------------------------------------------
data.list = readRDS(url("https://github.com/bblodfon/emba/blob/main/vignettes/data.rds?raw=true"))

model.predictions = data.list$model.predictions
models.stable.state = data.list$models.stable.state
models.link.operator = data.list$models.equations
observed.synergies = data.list$observed.synergies

# (x,y) coordinates for visualization
nice.layout = data.list$nice.layout
# model network as an igraph object
net = data.list$net

# drug combinations
drug.combos = colnames(model.predictions)

# change model names (shorter names for readability)
model.names = paste0("model", 1:7500)
rownames(model.predictions) = model.names
rownames(models.stable.state) = model.names
rownames(models.link.operator) = model.names

## ----input-2------------------------------------------------------------------
model.predictions[1:5, 77:84] %>% kable(caption = "Model predictions example")

## ----input-3------------------------------------------------------------------
models.stable.state[1:5, 5:11] %>% kable(caption = "Model stable states example")

## ----input-4------------------------------------------------------------------
models.link.operator[1:5, 1:10] %>% kable(caption = "Models link operator example")

## ----input-5, results='asis'--------------------------------------------------
usefun::pretty_print_vector_values(observed.synergies, vector.values.str = "observed synergies")

## ----tp-analysis-1------------------------------------------------------------
tp.analysis.res = emba::biomarker_tp_analysis(
  model.predictions, 
  models.stable.state, 
  models.link.operator, 
  observed.synergies, 
  penalty = 0.1,
  threshold = 0.55)

## ----tp-analysis-2, results = 'asis'------------------------------------------
usefun::pretty_print_vector_values(tp.analysis.res$predicted.synergies, vector.values.str = "predicted synergies")

## ----tp-analysis-3, fig.align='center', fig.width=7, fig.height=5.6-----------
pr(emba::make_barplot_on_models_stats(table(tp.analysis.res$models.synergies.tp), 
  title = "True Positive Synergy Predictions",
  xlab = "Number of maximum correctly predicted synergies",
  ylab = "Number of models"))

## ----tp-analysis-4------------------------------------------------------------
tp.analysis.res$diff.state.tp.mat %>% 
  as.data.frame() %>%
  select(c("AKT","PTEN","PSEN1","STAT3","CEBPA")) %>% # show only part of the matrix
  kable(caption = "Average Activity Difference Matrix")

## ----tp-analysis-5, results='asis'--------------------------------------------
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.active,
  vector.values.str = "active biomarkers")
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.inhibited,
  vector.values.str = "inhibited biomarkers")

## ----tp-analysis-6, fig.align='center', fig.width=7, fig.height=5.5-----------
pr(emba::plot_avg_state_diff_graph(net, tp.analysis.res$diff.state.tp.mat["(2,3)",], 
  layout = nice.layout, title = "Bad models (2 TP) vs Good models (3 TP)"))

## ----tp-analysis-7, results='asis'--------------------------------------------
tp.analysis.res.biased = emba::biomarker_tp_analysis(
  model.predictions, 
  models.stable.state, 
  models.link.operator, 
  observed.synergies, 
  penalty = 0,
  threshold = 0.7)

usefun::pretty_print_vector_values(tp.analysis.res.biased$biomarkers.tp.active,
  vector.values.str = "active biomarkers")

usefun::pretty_print_vector_values(tp.analysis.res.biased$biomarkers.tp.inhibited,
  vector.values.str = "inhibited biomarkers")

## ----tp-analysis-8------------------------------------------------------------
tp.analysis.res$diff.link.tp.mat %>% 
  as.data.frame() %>%
  select(c("AKT","PTEN","PSEN1","STAT3","CEBPA")) %>% # show only part of the matrix
  kable(caption = "Average Link Operator Difference Matrix")

## ----tp-analysis-9, results='asis'--------------------------------------------
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.or,
  vector.values.str = "'OR' biomarkers")
usefun::pretty_print_vector_values(tp.analysis.res$biomarkers.tp.and,
  vector.values.str = "'AND' biomarkers")

## ----tp-analysis-10, fig.align='center', fig.width=7, fig.height=5.5----------
pr(emba::plot_avg_link_operator_diff_graph(net, tp.analysis.res$diff.link.tp.mat["(2,3)",], 
  layout = nice.layout, title = "Bad models (2 TP) vs Good models (3 TP)"))

## ----mcc-analysis-1-----------------------------------------------------------
mcc.analysis.res = emba::biomarker_mcc_analysis(
  model.predictions, 
  models.stable.state, 
  models.link.operator, 
  observed.synergies, 
  threshold = 0.65,
  num.of.mcc.classes = 4,
  penalty = 0.2)

## ----mcc-analysis-2, results = 'asis'-----------------------------------------
usefun::pretty_print_vector_values(mcc.analysis.res$predicted.synergies, vector.values.str = "predicted synergies")

## ----mcc-analysis-3, fig.align='center', fig.width=7, fig.height=5.8----------
pr(emba::make_barplot_on_models_stats(table(mcc.analysis.res$models.mcc), 
  title = "MCC scores", xlab = "MCC value", 
  ylab = "Number of models", cont.values = TRUE))

## ----mcc-analysis-4, fig.align='center', fig.width=7, fig.height=5.5----------
models.mcc = mcc.analysis.res$models.mcc
num.of.mcc.classes = 4

res = Ckmeans.1d.dp(x = models.mcc, k = num.of.mcc.classes)
models.cluster.ids = res$cluster

pr(emba::plot_mcc_classes_hist(models.mcc, models.cluster.ids, num.of.mcc.classes))

## ----mcc-analysis-5-----------------------------------------------------------
mcc.analysis.res$diff.state.mcc.mat %>%
  as.data.frame() %>%
  select(c("AKT","PPM1A","PTEN","PSEN1","PTK2","CEBPA")) %>% # show only part of the matrix
  kable(caption = "Average Activity Difference Matrix")

## ----mcc-analysis-6, results='asis'-------------------------------------------
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.active,
  vector.values.str = "active biomarkers")
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.inhibited,
  vector.values.str = "inhibited biomarkers")

## ----mcc-analysis-7, fig.align='center', fig.width=7, fig.height=5.5----------
pr(emba::plot_avg_state_diff_graph(net, mcc.analysis.res$diff.state.mcc.mat["(1,4)",], 
  layout = nice.layout, title = "Bad models (MCC Class 1) vs Good models (MCC Class 4)"))

## ----mcc-analysis-8-----------------------------------------------------------
mcc.analysis.res$diff.link.mcc.mat %>% 
  as.data.frame() %>%
  select(c("AKT","PTEN","PSEN1","CEBPA","STAT3","JAK1")) %>% # show only part of the matrix
  kable(caption = "Average Link Operator Difference Matrix")

## ----mcc-analysis-9, results='asis'-------------------------------------------
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.or,
  vector.values.str = "'OR' biomarkers")
usefun::pretty_print_vector_values(mcc.analysis.res$biomarkers.mcc.and,
  vector.values.str = "'AND' biomarkers")

## ----mcc-analysis-10, fig.align='center', fig.width=7, fig.height=5.5---------
pr(emba::plot_avg_link_operator_diff_graph(net, mcc.analysis.res$diff.link.mcc.mat["(1,4)",], 
  layout = nice.layout, title = "Bad models (MCC Class 1) vs Good models (MCC Class 4)"))

## ----synergy-analysis-1-------------------------------------------------------
synergy.analysis.res = emba::biomarker_synergy_analysis(
  model.predictions,
  models.stable.state,
  models.link.operator,
  observed.synergies,
  threshold = 0.5,
  calculate.subsets.stats = TRUE,
  penalty = 0.1)

## ----synergy-analysis-2, fig.align='center', fig.width=7, fig.height=5.5------
pr(emba::make_barplot_on_synergy_subset_stats(
  synergy.analysis.res$synergy.subset.stats,
  threshold.for.subset.removal = 1, 
  bottom.mar = 9))

## ----synergy-analysis-3, fig.align='center', fig.width=7, fig.height=5--------
synergy.analysis.res$diff.state.synergies.mat[,1:8] %>% 
  kable(caption = "Average State Differences per Synergy Predicted", digits = 3)

## ----synergy-analysis-4, fig.align='center', fig.width=7, fig.height=5.5------
pr(emba::plot_avg_state_diff_graph(net,
  synergy.analysis.res$diff.state.synergies.mat["PI-D1",],
  layout = nice.layout, title = "Prediction of PI-D1 (Good Models: synergy, Bad Models: antagonism)"))

## ----synergy-analysis-5, fig.align='center', fig.width=7, fig.height=5--------
# prune nodes (columns) that were not found as biomarkers for any predicted synergy
biomarker.act.mat = usefun::prune_columns_from_df(
  df = synergy.analysis.res$activity.biomarkers, value = 0)

biomarker.act.mat[, 4:12] %>% # show only part of the matrix
  kable(caption = "Activity State Biomarkers Per Synergy Predicted")

## ---- eval = FALSE------------------------------------------------------------
#  # define your own threshold
#  my.thres = 0.76
#  activity.biomarkers.new = as.data.frame(apply(
#    synergy.analysis.res$diff.state.synergies.mat, c(1,2),
#    usefun::get_ternary_class_id, my.thres))

## ----synergy-analysis-6, results='asis'---------------------------------------
drug.comb = "AK-PD"

syn.models.num = sum(model.predictions[, drug.comb] == 1 & !is.na(model.predictions[, drug.comb]))
ant.models.num  = sum(model.predictions[, drug.comb] == 0 & !is.na(model.predictions[, drug.comb]))

usefun::pretty_print_string(paste0("Number of models (AK-PD): #Synergistic: ", syn.models.num, ", #Antagonistic: ", ant.models.num))

## ----r-session-info, comment=""-----------------------------------------------
sessionInfo()

