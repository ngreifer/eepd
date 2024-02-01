# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
library(tidyverse)
library(ggthemes)
library(magrittr)
library(formula.tools)
load(file = "_dev/data.RData")


# Model Selection ---------------------------------------------------------

#family <- c("none", "poisson", "quasipoisson", "negative binomial")
outcome_diff <- c("none", "first")
outcome_log <- c("No", "Yes")
time_trend <- c("None", "Linear", "Quadratic", "Cubic")
lagged_DV <- c("No", "Yes")
mod_combs <- expand.grid(outcome_diff, outcome_log, time_trend, lagged_DV)
colnames(mod_combs) <- c("outcome_diff", "outcome_log", "time_trend", "lagged_DV")
mod_combs

val_years <- 1999:2004
val_index <- sort(x = unique(data$time_index[which(data$year %in% val_years)]),
                  decreasing = TRUE)

source("eepd_fun_sim.R")

## Exclude models with both time_index and time FEs
model_args <- list(list("formula" = crude_rate ~ 1, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + time_index, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + time_index, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + time_index, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + time_index, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + time_index, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + time_index, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + time_index_squared, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + time_index_squared, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + time_index_squared, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + time_index_squared, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + time_index_squared, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + time_index_squared, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + time_index_cubed, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + time_index_cubed, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + time_index_cubed, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + time_index_cubed, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + time_index_cubed, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + time_index_cubed, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + crude_rate_lag_1, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + crude_rate_lag_1, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + crude_rate_lag_1, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + crude_rate_lag_1, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + time_index + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + time_index + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + time_index + crude_rate_lag_1, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + time_index + crude_rate_lag_1, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + time_index + crude_rate_lag_1, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + time_index + crude_rate_lag_1, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + time_index_squared + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + time_index_squared + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + time_index_squared + crude_rate_lag_1, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + time_index_squared + crude_rate_lag_1, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + time_index_squared + crude_rate_lag_1, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + time_index_squared + crude_rate_lag_1, "family" = NULL, diff_k = 1, log = FALSE),
                   
                   list("formula" = crude_rate ~ 1 + time_index_cubed + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = FALSE),
                   list("formula" = crude_rate ~ 1 + time_index_cubed + crude_rate_lag_1, "family" = NULL, diff_k = 0, log = TRUE),
                   #list("formula" = crude_rate ~ 1 + time_index_cubed + crude_rate_lag_1, "family" = "poisson", diff_k = 0, log = FALSE),
                   #list("formula" = crude_rate ~ 1 + time_index_cubed + crude_rate_lag_1, "family" = "quasipoisson", diff_k = 0),
                   #list("formula" = crude_rate ~ 1 + time_index_cubed + crude_rate_lag_1, "family" = "negative binomial", diff_k = 0),
                   
                   list("formula" = crude_rate ~ 1 + time_index_cubed + crude_rate_lag_1, "family" = NULL, diff_k = 1, log = FALSE))

n_sims <- 5

model_sel_out <- matrix(data = NA,
                        nrow = length(val_index),
                        ncol = length(model_args),
                        dimnames = list(paste("val_time", val_index, sep = "_"),
                                        paste("model", 1:length(model_args), sep = "_")))

model_sel_out_list <- lapply(X = seq_len(n_sims),
                             FUN = function(X) { model_sel_out })
library(MASS)
set.seed(11242017)

for(j in 1:n_sims){
  cat(j, " ")
  for(i in val_index){
    
    model_sel_out_list[[j]][paste("val_time", i, sep = "_"),] = sapply(X = 1:length(model_args),
                                                                       FUN = function(x) { 
                                                                         eepd_sim(model = model_args[[x]]$formula,
                                                                                  family = model_args[[x]]$family,
                                                                                  group_var = group,
                                                                                  unit_var = state,
                                                                                  time_var = time,
                                                                                  data = data,
                                                                                  val_time = i,
                                                                                  diff_k = model_args[[x]]$diff_k,
                                                                                  log = model_args[[x]]$log)$avg_diff_pred_errors })
  }
  
}

## min max loss function
max_losses <- lapply(X = 1:length(model_sel_out_list),
                     FUN = function(j) { apply(X = model_sel_out_list[[j]],
                                               MARGIN = 2,
                                               FUN = function(x) { max(abs(x)) }) })

# ATT Estimation ----------------------------------------------------------
source("eepd_fun_est_ATT.R")
est_ATTs <- sapply(X = 1:length(max_losses),
                   FUN = function(x) { eepd_est_ATT(model = model_args[[which.min(max_losses[[x]])]]$formula,
                                                    family = model_args[[which.min(max_losses[[x]])]]$family,
                                                    group_var = group,
                                                    unit_var = state,
                                                    time_var = time,
                                                    data = data,
                                                    val_time = 15,
                                                    diff_k = model_args[[which.min(max_losses[[x]])]]$diff_k,
                                                    log = model_args[[which.min(max_losses[[x]])]]$log)$avg_diff_pred_errors })
save(est_ATTs, file = "est_ATTs.RData")



