# Setup -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
library(tidyverse)
library(ggthemes)
library(magrittr)
library(pbapply)
library(parallel)
library(MASS)
library(formula.tools)
data_92_98 <- read.delim(file = "_dev/all_states_pre_period.txt",
                         header = TRUE,
                         sep = "\t",
                         dec = ".")
data_99_07 <- read.delim(file = "_dev/all_states_before_after_periods.txt",
                         header = TRUE,
                         sep = "\t",
                         dec = ".") 

# Data Cleaning -----------------------------------------------------------

## Hasewaga et al (2019), p. 375
treated_state <- "Missouri"
control_states <- c("Arkansas",
                    "Illinois",
                    "Iowa",
                    "Kansas",
                    "Kentucky",
                    "Nebraska",
                    "Oklahoma",
                    "Tennessee")
states <- union(x = treated_state,
                y = control_states)
## Training data: 1994 - 1998
train_years <- 1994:1998
## Validation data: 1999 - 2007
val_years <- 1999:2007
## Post-treatment data: 2008 - 2016
post_treat_years <- 2008:2016

data <- rbind(data_92_98, data_99_07)
## Use only treated state and 8 neighboring control states, as well as years after 1994
data <- dplyr::filter(.data = data, State %in% states & Year >= 1994) %>%
    dplyr::select(State, Year, Deaths, Crude.Rate, Age.Adjusted.Rate) %>%
    dplyr::rename(state = State,
                  year = Year,
                  deaths = Deaths,
                  crude_rate = Crude.Rate,
                  age_adj_rate = Age.Adjusted.Rate)
data <- mutate(.data = data,
               group = ifelse(test = state == "Missouri",
                              yes = 1,
                              no = 0),
               time = year - (min(year) - 1),
               time_index = year - (min(year) - 1),
               time_index_squared = time_index^2,
               time_index_cubed = time_index^3,
               group = ifelse(test = state == "Missouri",
                              yes = 1,
                              no = 0),
               crude_rate = gsub(pattern = "[^0-9.-]", ## Clean up outcomes, two of which described as "Unreliable"
                                 replacement = "",
                                 x = crude_rate),
               crude_rate = as.numeric(crude_rate),
               log_crude_rate = log(crude_rate),
               age_adj_rate = gsub(pattern = "[^0-9.-]", ## Clean up outcomes, two of which described as "Unreliable"
                                   replacement = "",
                                   x = age_adj_rate),
               age_adj_rate = as.numeric(age_adj_rate),
               log_age_adj_rate = log(age_adj_rate),
               state = as.factor(state),
               treat = ifelse(test = state == "Missouri" & year %in% post_treat_years,
                              yes = 1,
                              no = 0)) %>%
    arrange(state, year, group)

data <- group_by(.data = data,
                 state) %>%
    mutate(crude_rate_lag_1 = lag(x = crude_rate,
                                  n = 1,
                                  default = NA),
           crude_rate_lag_2 = lag(x = crude_rate,
                                  n = 2,
                                  default = NA),
           age_adj_rate_lag_1 = lag(x = age_adj_rate,
                                    n = 1,
                                    default = NA),
           age_adj_rate_lag_2 = lag(x = age_adj_rate,
                                    n = 2,
                                    default = NA),
           time_index_squared = time_index^2,
           time_index_cubes = time_index^3)

## Recode first post-treatment period (time_index = 15) to be average outcome over all post-treatment periods
treated_post_treat_means <- filter(.data = data,
                                   group == 1) %>%
    group_by(state) %>%
    summarize(post_treat_mean = mean(crude_rate[which(year %in% post_treat_years)]))

control_post_treat_means <- filter(.data = data,
                                   group == 0) %>%
    group_by(state) %>%
    summarize(post_treat_mean = mean(crude_rate[which(year %in% post_treat_years)]))

data <- mutate(.data = data,
               crude_rate = ifelse(test = state == "Arkansas" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Arkansas"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Illinois" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Illinois"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Iowa" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Iowa"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Kansas" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Kansas"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Kentucky" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Kentucky"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Missouri" & time_index == 15,
                                   yes = treated_post_treat_means$post_treat_mean[treated_post_treat_means$state == "Missouri"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Nebraska" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Nebraska"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Oklahoma" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Oklahoma"],
                                   no = crude_rate),
               crude_rate = ifelse(test = state == "Tennessee" & time_index == 15,
                                   yes = control_post_treat_means$post_treat_mean[control_post_treat_means$state == "Tennessee"],
                                   no = crude_rate)) %>%
    filter(time_index <= 15)

# Estimation ---------------------------------------------------------------

models <- list("model_a" = (crude_rate ~ 1),
               "model_b" = (crude_rate ~ 1 + time_index))

pred_cor_indiv <- function(model,
                           unit_var,
                           time_var,
                           data,
                           unit,
                           val_time,
                           sim = FALSE){
    
    indiv_train_data = filter(.data = data,
                              {{time_var}} < val_time & {{unit_var}} == unit)
    indiv_val_data = filter(.data = data,
                            {{time_var}} == val_time & {{unit_var}} == unit)
    
    mod = lm(formula = model,
             data = indiv_train_data)
    
    mod_mat = model.matrix(lm(formula = model,
                              data = indiv_val_data))
    
    indiv_outcome = pull(.data = indiv_val_data, var = formula.tools::get.vars(model)[1])
    
    return(ifelse(sim,
                  yes = indiv_outcome - as.numeric(mod_mat %*% MASS::mvrnorm(n = 1, mu = coef(mod), Sigma = vcov(mod))),
                  no = indiv_outcome - as.numeric(mod_mat %*% coef(mod))))
    
}
ATT_est_model_a <- pred_cor_indiv(model = models[["model_a"]],
                                  unit_var = state,
                                  time_var = time,
                                  data = data,
                                  unit = "Missouri",
                                  val_time = 15,
                                  sim = FALSE) - mean(pbsapply(X = control_states,
                                                               FUN = function(x) { pred_cor_indiv(model = models[["model_a"]],
                                                                                                  unit_var = state,
                                                                                                  time_var = time,
                                                                                                  data = data,
                                                                                                  unit = x,
                                                                                                  val_time = 15,
                                                                                                  sim = FALSE) },
                                                               cl = NULL))
ATT_est_model_b <- pred_cor_indiv(model = models[["model_b"]],
                                  unit_var = state,
                                  time_var = time,
                                  data = data,
                                  unit = "Missouri",
                                  val_time = 15,
                                  sim = FALSE) - mean(pbsapply(X = control_states,
                                                               FUN = function(x) { pred_cor_indiv(model = models[["model_b"]],
                                                                                                  unit_var = state,
                                                                                                  time_var = time,
                                                                                                  data = data,
                                                                                                  unit = x,
                                                                                                  val_time = 15,
                                                                                                  sim = FALSE) },
                                                               cl = NULL))

# Laura's approach ----------------------------------------

n_boots <- 10

boot_ATT_ests_model_a <- rep(x = NA, times = n_boots)
boot_ATT_ests_model_b <- rep(x = NA, times = n_boots)

set.seed(11242017)

for(i in 1:n_boots){
    
    boot_control_states = sample(x = control_states, replace = TRUE)
    
    boot_ATT_ests_model_a[i] = pred_cor_indiv(model = models[["model_a"]],
                                              unit_var = state,
                                              time_var = time,
                                              data = data,
                                              unit = "Missouri",
                                              val_time = 15,
                                              sim = FALSE) - mean(pbsapply(X = boot_control_states,
                                                                           FUN = function(x) { pred_cor_indiv(model = models[["model_a"]],
                                                                                                              unit_var = state,
                                                                                                              time_var = time,
                                                                                                              data = data,
                                                                                                              unit = x,
                                                                                                              val_time = 15,
                                                                                                              sim = FALSE) },
                                                                           cl = NULL))
    boot_ATT_ests_model_b[i] = pred_cor_indiv(model = models[["model_b"]],
                                              unit_var = state,
                                              time_var = time,
                                              data = data,
                                              unit = "Missouri",
                                              val_time = 15,
                                              sim = FALSE) - mean(pbsapply(X = boot_control_states,
                                                                           FUN = function(x) { pred_cor_indiv(model = models[["model_b"]],
                                                                                                              unit_var = state,
                                                                                                              time_var = time,
                                                                                                              data = data,
                                                                                                              unit = x,
                                                                                                              val_time = 15,
                                                                                                              sim = FALSE) },
                                                                           cl = NULL))
}

val_index <- sort(x = unique(data$time_index[which(data$year %in% val_years)]),
                  decreasing = TRUE)

n_sims <- 10
min_max_out <- matrix(data = NA,
                      nrow = length(models),
                      ncol = n_sims,
                      dimnames = list(names(models)))

set.seed(11242017)
for(i in 1:n_sims){
    
    min_max_out["model_a", i] = max(abs(pbsapply(X = val_index,
                                                 FUN = function(t){
                                                     
                                                     pred_cor_indiv(model = models[["model_a"]],
                                                                    unit_var = state,
                                                                    time_var = time,
                                                                    data = data,
                                                                    unit = "Missouri",
                                                                    val_time = t,
                                                                    sim = TRUE) - mean(pbsapply(X = control_states,
                                                                                                FUN = function(x) { pred_cor_indiv(model = models[["model_a"]],
                                                                                                                                   unit_var = state,
                                                                                                                                   time_var = time,
                                                                                                                                   data = data,
                                                                                                                                   unit = x,
                                                                                                                                   val_time = t,
                                                                                                                                   sim = TRUE) },
                                                                                                cl = NULL)) },
                                                 cl = NULL)))
    
    min_max_out["model_b", i] = max(abs(pbsapply(X = val_index,
                                                 FUN = function(t){
                                                     
                                                     pred_cor_indiv(model = models[["model_b"]],
                                                                    unit_var = state,
                                                                    time_var = time,
                                                                    data = data,
                                                                    unit = "Missouri",
                                                                    val_time = t,
                                                                    sim = TRUE) - mean(pbsapply(X = control_states,
                                                                                                FUN = function(x) { pred_cor_indiv(model = models[["model_b"]],
                                                                                                                                   unit_var = state,
                                                                                                                                   time_var = time,
                                                                                                                                   data = data,
                                                                                                                                   unit = x,
                                                                                                                                   val_time = t,
                                                                                                                                   sim = TRUE) },
                                                                                                cl = NULL)) },
                                                 cl = NULL)))
}

opt_mods_out <- apply(X = min_max_out,
                      MARGIN = 2,
                      FUN = function(x) { row.names(min_max_out)[which.min(x)] })
opt_mod_probs <- sapply(X = row.names(min_max_out),
                        FUN = function(x) { mean(opt_mods_out == x) })

weighted_boot_ests <- boot_ATT_ests_model_a * opt_mod_probs[1] + boot_ATT_ests_model_b * opt_mod_probs[2]
var(weighted_boot_ests)
quantile(x = weighted_boot_ests, probs = c(0.025, 0.975))

# New approach ------------------------------------------------------------

n_boots <- 10
boot_ATT_ests <- vector(mode = "list", length = n_boots)

n_sims <- 10

set.seed(11242017)
for(i in 1:n_boots){
    
    boot_control_states = sample(x = control_states, replace = TRUE)
    
    for(j in 1:n_sims){
        
        max_abs_errs = c(max(abs(pbsapply(X = val_index,
                                          FUN = function(t){
                                              
                                              pred_cor_indiv(model = models[["model_a"]],
                                                             unit_var = state,
                                                             time_var = time,
                                                             data = data,
                                                             unit = "Missouri",
                                                             val_time = t,
                                                             sim = TRUE) -
                                                  mean(pbsapply(X = boot_control_states,
                                                                function(x) {
                                                                    pred_cor_indiv(model = models[["model_a"]],
                                                                                   unit_var = state,
                                                                                   time_var = time,
                                                                                   data = data,
                                                                                   unit = x,
                                                                                   val_time = t,
                                                                                   sim = TRUE)
                                                                },
                                                                cl = NULL)) },
                                          cl = NULL))),
                         max(abs(pbsapply(X = val_index,
                                          FUN = function(t){
                                              
                                              pred_cor_indiv(model = models[["model_b"]],
                                                             unit_var = state,
                                                             time_var = time,
                                                             data = data,
                                                             unit = "Missouri",
                                                             val_time = t,
                                                             sim = TRUE) - mean(pbsapply(X = boot_control_states,
                                                                                         FUN = function(x) { pred_cor_indiv(model = models[["model_b"]],
                                                                                                                            unit_var = state,
                                                                                                                            time_var = time,
                                                                                                                            data = data,
                                                                                                                            unit = x,
                                                                                                                            val_time = t,
                                                                                                                            sim = TRUE) },
                                                                                         cl = NULL)) },
                                          cl = NULL))))
        
        boot_ATT_ests[[i]][j] = pred_cor_indiv(model = models[[names(models)[which.min(max_abs_errs)]]],
                                               unit_var = state,
                                               time_var = time,
                                               data = data,
                                               unit = "Missouri",
                                               val_time = 15,
                                               sim = FALSE) - mean(pbsapply(X = boot_control_states,
                                                                            FUN = function(x) { pred_cor_indiv(model = models[[names(models)[which.min(max_abs_errs)]]],
                                                                                                               unit_var = state,
                                                                                                               time_var = time,
                                                                                                               data = data,
                                                                                                               unit = x,
                                                                                                               val_time = 15,
                                                                                                               sim = FALSE) },
                                                                            cl = NULL))
        
    }
    
}

var(unlist(boot_ATT_ests))
quantile(x = unlist(boot_ATT_ests), probs = c(0.025, 0.975))
