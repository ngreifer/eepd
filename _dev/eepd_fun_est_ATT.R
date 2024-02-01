#' Equal-expected-prediction-errors placebo test
#' This function returns the absolute difference in average prediction errors between treated and control groups
#' @param model a formula object, the same type that would be fed to formula argument of lm()
#' @param group_var the variable in dataset (entered without quotation marks) that denotes the binary (treatment vs control) group
#' @param unit_var the variable in dataset (entered without quotation marks) for unit fixed effects
#' @param time_var the variable in dataset (entered without quotation marks) for time fixed effects
#' @param unit_FEs TRUE or FALSE indicating whether to use unit fixed effects or not
#' @param time_FEs TRUE or FALSE indicating whether to use time fixed effects or not
#' @param data the dataset where the unit_var and time_var variables are located, as well as variables in model argument, akin to data argument of lm()
#' @param val_time the integer-valued time period (of the time_var variable) that serves as the period for validation data; training data is all periods before val_time
eepd_est_ATT <- function(model,
                         family = NULL, ## "poisson", "quasipoisson" or "negative binomial"
                         group_var,
                         unit_var,
                         time_var,
                         data,
                         val_time,
                         diff_k = 0, ## k is the number used for differencing, 0 is original, 1 is first diff, etc
                         log = FALSE){ 
  
  ## Log of outcome on original scale
  data = mutate(.data = data,
                dplyr::across(.cols = formula.tools::get.vars(model)[1],
                              .fns = list(log = ~ log(.x))))
  
  ## Create lagged variables on original scale
  if(diff_k > 0){
    
    data = mutate(.data = data,
                  dplyr::across(.cols = formula.tools::get.vars(model)[1],
                                .fns = list(~ .x - lag(x = .x,
                                                       n = diff_k,
                                                       default = NA)),
                                .names = paste(formula.tools::get.vars(model)[1],
                                               "diff",
                                               diff_k,
                                               sep = "_")))
    
    ## Create lagged variables on log scale
    data = mutate(.data = data,
                  dplyr::across(.cols = paste(formula.tools::get.vars(model)[1], "log", sep = "_"),
                                .fns = list(~ .x - lag(x = .x,
                                                       n = diff_k,
                                                       default = NA)),
                                .names = paste(paste(formula.tools::get.vars(model)[1],
                                                     "log",
                                                     sep = "_"),
                                               "diff",
                                               diff_k,
                                               sep = "_")))
    
    
  } else{
    
    data = mutate(.data = data,
                  dplyr::across(.cols = formula.tools::get.vars(model)[1],
                                .fns = list(~ .x ),
                                .names = paste(formula.tools::get.vars(model)[1],
                                               "diff",
                                               diff_k,
                                               sep = "_")))
    
    data = mutate(.data = data,
                  dplyr::across(.cols = paste(formula.tools::get.vars(model)[1], "log", sep = "_"),
                                .fns = list(~ .x ),
                                .names = paste(paste(formula.tools::get.vars(model)[1],
                                                     "log",
                                                     sep = "_"),
                                               "diff",
                                               diff_k,
                                               sep = "_")))
    
  }
  
  ## Select training data (all periods less than the validation period)
  train_data = filter(.data = data,
                      {{time_var}} < val_time) ## {{col}} is a string literal where expression in {{}} is evaluated and then inserted into argument string
  
  ## Select validation data
  val_data = filter(.data = data,
                    {{time_var}} <= val_time)
  
  treated_train_data = filter(.data = train_data,
                              {{time_var}} < val_time & {{group_var}} == 1)
  treated_val_data = filter(.data = val_data,
                            {{time_var}} == val_time & {{group_var}} == 1)
  
  control_train_data = filter(.data = train_data,
                              {{time_var}} < val_time & {{group_var}} == 0)
  control_val_data = filter(.data = val_data,
                            {{time_var}} == val_time & {{group_var}} == 0)
  
  if(log == FALSE){
    
    if(length(formula.tools::get.vars(model)) > 1) { new_form_names = c(paste(formula.tools::get.vars(model)[1], "diff", diff_k, sep = "_"),
                                                                        formula.tools::get.vars(model)[-1]) } else{
                                                                          
                                                                          new_form_names = paste(formula.tools::get.vars(model)[1], "diff", diff_k, sep = "_")
                                                                        }
    
    ## If model has no predictors (just an intercept) then argument within if() statement will not work, so use if else
    if(length(new_form_names) > 1) { new_form = as.formula(paste(new_form_names[1], "~", paste(new_form_names[-1], collapse = "+")))} else{
      new_form = as.formula(paste(new_form_names[1], "~", "1")) }
    
    if(is.null(family)){
      
      treated_train_fits = treated_train_data %>%
        group_by({{unit_var}}) %>%
        do(mod = lm(formula = new_form,
                    data = .))
      control_train_fits = control_train_data %>%
        group_by({{unit_var}}) %>%
        do(mod = lm(formula = new_form,
                    data = .))
      
      treated_val_mod_mats = treated_val_data %>%
        group_by({{unit_var}}) %>%
        do(mod_mat = model.matrix(lm(formula = new_form,
                                     data = .)))
      control_val_mod_mats = control_val_data %>%
        group_by({{unit_var}}) %>%
        do(mod_mat = model.matrix(lm(formula = new_form,
                                     data = .)))
      
      treated_val_outcomes = pull(.data = treated_val_data, var = formula.tools::get.vars(model)[1])
      control_val_outcomes = pull(.data = control_val_data, var = formula.tools::get.vars(model)[1])
      
      treated_val_preds_trans_scale = sapply(X = 1:nrow(treated_train_fits),
                                             FUN = function(x) { 
                                               as.numeric(treated_val_mod_mats[["mod_mat"]][[x]] %*% coef(treated_train_fits[["mod"]][[x]])) })
      control_val_preds_trans_scale = sapply(X = 1:nrow(control_train_fits), 
                                             FUN = function(x) { 
                                               as.numeric(control_val_mod_mats[["mod_mat"]][[x]] %*% coef(control_train_fits[["mod"]][[x]])) })
      
      if(diff_k == 0){ 
        
        treated_val_preds = treated_val_preds_trans_scale
        control_val_preds = control_val_preds_trans_scale
        
      } else{
        
        treated_val_preds = (treated_val_preds_trans_scale) +
          pull(.data = filter(.data = treated_train_data, {{time_var}} == (val_time - diff_k)),
               var = formula.tools::get.vars(model)[1])
        control_val_preds = (control_val_preds_trans_scale) +
          pull(.data = filter(.data = control_train_data, {{time_var}} == (val_time - diff_k)),
               var = formula.tools::get.vars(model)[1])
      }
      
    } else{
      
      if(family == "poisson"){
        
        treated_train_fits = treated_train_data %>%
          group_by({{unit_var}}) %>%
          do(mod = glm(formula = new_form,
                       family = poisson(link = "log"),
                       data = .))
        control_train_fits = control_train_data %>%
          group_by({{unit_var}}) %>%
          do(mod = glm(formula = new_form,
                       family = poisson(link = "log"),
                       data = .))
        
        treated_val_mod_mats = treated_val_data %>%
          group_by({{unit_var}}) %>%
          do(mod_mat = model.matrix(glm(formula = new_form,
                                        family = poisson(link = "log"),
                                        data = .)))
        control_val_mod_mats = control_val_data %>%
          group_by({{unit_var}}) %>%
          do(mod_mat = model.matrix(glm(formula = new_form,
                                        family = poisson(link = "log"),
                                        data = .)))
        
        treated_val_outcomes = pull(.data = treated_val_data, var = formula.tools::get.vars(model)[1])
        control_val_outcomes = pull(.data = control_val_data, var = formula.tools::get.vars(model)[1])
        
        treated_val_preds_trans_scale = sapply(X = 1:nrow(treated_train_fits),
                                               FUN = function(x) { 
                                                 as.numeric(treated_val_mod_mats[["mod_mat"]][[x]] %*% coef(treated_train_fits[["mod"]][[x]])) })
        control_val_preds_trans_scale = sapply(X = 1:nrow(control_train_fits), 
                                               FUN = function(x) { 
                                                 as.numeric(control_val_mod_mats[["mod_mat"]][[x]] %*% coef(control_train_fits[["mod"]][[x]])) })
        
        
        if(diff_k == 0){ 
          
          treated_val_preds = treated_val_preds_trans_scale
          control_val_preds = control_val_preds_trans_scale
          
        } else{
          
          treated_val_preds = (treated_val_preds_trans_scale) +
            pull(.data = filter(.data = treated_train_data, {{time_var}} == (val_time - diff_k)),
                 var = formula.tools::get.vars(model)[1])
          control_val_preds = (control_val_preds_trans_scale) +
            pull(.data = filter(.data = control_train_data, {{time_var}} == (val_time - diff_k)),
                 var = formula.tools::get.vars(model)[1])
        }
        
      } else{
        
        if(family == "quasipoisson"){
          
          treated_train_fits = treated_train_data %>%
            group_by({{unit_var}}) %>%
            do(mod = glm(formula = new_form,
                         family = quasipoisson(link = "log"),
                         data = .))
          control_train_fits = control_train_data %>%
            group_by({{unit_var}}) %>%
            do(mod = glm(formula = new_form,
                         family = quasipoisson(link = "log"),
                         data = .))
          
          treated_val_mod_mats = treated_val_data %>%
            group_by({{unit_var}}) %>%
            do(mod_mat = model.matrix(glm(formula = new_form,
                                          family = quasipoisson(link = "log"),
                                          data = .)))
          control_val_mod_mats = control_val_data %>%
            group_by({{unit_var}}) %>%
            do(mod_mat = model.matrix(glm(formula = new_form,
                                          family = quasipoisson(link = "log"),
                                          data = .)))
          
          treated_val_outcomes = pull(.data = treated_val_data, var = formula.tools::get.vars(model)[1])
          control_val_outcomes = pull(.data = control_val_data, var = formula.tools::get.vars(model)[1])
          
          treated_val_preds_trans_scale = sapply(X = 1:nrow(treated_train_fits),
                                                 FUN = function(x) { 
                                                   as.numeric(treated_val_mod_mats[["mod_mat"]][[x]] %*% coef(treated_train_fits[["mod"]][[x]])) })
          control_val_preds_trans_scale = sapply(X = 1:nrow(control_train_fits), 
                                                 FUN = function(x) { 
                                                   as.numeric(control_val_mod_mats[["mod_mat"]][[x]] %*% coef(control_train_fits[["mod"]][[x]])) })
          
          
          if(diff_k == 0){ 
            
            treated_val_preds = treated_val_preds_trans_scale
            control_val_preds = control_val_preds_trans_scale
            
          } else{
            
            treated_val_preds = (treated_val_preds_trans_scale) +
              pull(.data = filter(.data = treated_train_data, {{time_var}} == (val_time - diff_k)),
                   var = formula.tools::get.vars(model)[1])
            control_val_preds = (control_val_preds_trans_scale) +
              pull(.data = filter(.data = control_train_data, {{time_var}} == (val_time - diff_k)),
                   var = formula.tools::get.vars(model)[1])
          }
          
        } else{ ## negative binomial model
          
          treated_train_fits = treated_train_data %>%
            group_by({{unit_var}}) %>%
            do(mod = MASS::glm.nb(formula = new_form,
                                  data = .))
          control_train_fits = control_train_data %>%
            group_by({{unit_var}}) %>%
            do(mod = MASS::glm.nb(formula = new_form,
                                  data = .))
          
          treated_val_mod_mats = treated_val_data %>%
            group_by({{unit_var}}) %>%
            do(mod_mat = model.matrix(MASS::glm.nb(formula = new_form,
                                                   data = .)))
          control_val_mod_mats = control_val_data %>%
            group_by({{unit_var}}) %>%
            do(mod_mat = model.matrix(MASS::glm.nb(formula = new_form,
                                                   data = .)))
          
          treated_val_outcomes = pull(.data = treated_val_data, var = formula.tools::get.vars(model)[1])
          control_val_outcomes = pull(.data = control_val_data, var = formula.tools::get.vars(model)[1])
          
          treated_val_preds_trans_scale = sapply(X = 1:nrow(treated_train_fits),
                                                 FUN = function(x) { 
                                                   as.numeric(treated_val_mod_mats[["mod_mat"]][[x]] %*% coef(treated_train_fits[["mod"]][[x]])) })
          control_val_preds_trans_scale = sapply(X = 1:nrow(control_train_fits), 
                                                 FUN = function(x) { 
                                                   as.numeric(control_val_mod_mats[["mod_mat"]][[x]] %*% coef(control_train_fits[["mod"]][[x]])) })
          
          if(diff_k == 0){ 
            
            treated_val_preds = treated_val_preds_trans_scale
            control_val_preds = control_val_preds_trans_scale
            
          } else{
            
            treated_val_preds = (treated_val_preds_trans_scale) +
              pull(.data = filter(.data = treated_train_data, {{time_var}} == (val_time - diff_k)),
                   var = formula.tools::get.vars(model)[1])
            control_val_preds = (control_val_preds_trans_scale) +
              pull(.data = filter(.data = control_train_data, {{time_var}} == (val_time - diff_k)),
                   var = formula.tools::get.vars(model)[1])
          }
          
        }
        
      }
      
      
    } } else{
      
      
      if(length(formula.tools::get.vars(model)) > 1) { new_form_names = c(paste(formula.tools::get.vars(model)[1], "log", "diff", diff_k, sep = "_"),
                                                                          formula.tools::get.vars(model)[-1]) } else{
                                                                            
                                                                            new_form_names = paste(formula.tools::get.vars(model)[1], "log", "diff", diff_k, sep = "_")
                                                                          }
      
      ## If model has no predictors (just an intercept) then argument within if() statement will not work, so use if else
      if(length(new_form_names) > 1) { new_form = as.formula(paste(new_form_names[1], "~", paste(new_form_names[-1], collapse = "+")))} else{
        new_form = as.formula(paste(new_form_names[1], "~", "1")) }
      
      treated_train_fits = treated_train_data %>%
        group_by({{unit_var}}) %>%
        do(mod = lm(formula = new_form,
                    data = .))
      control_train_fits = control_train_data %>%
        group_by({{unit_var}}) %>%
        do(mod = lm(formula = new_form,
                    data = .))
      
      treated_val_mod_mats = treated_val_data %>%
        group_by({{unit_var}}) %>%
        do(mod_mat = model.matrix(lm(formula = new_form,
                                     data = .)))
      control_val_mod_mats = control_val_data %>%
        group_by({{unit_var}}) %>%
        do(mod_mat = model.matrix(lm(formula = new_form,
                                     data = .)))
      
      treated_val_outcomes = pull(.data = treated_val_data, var = formula.tools::get.vars(model)[1])
      control_val_outcomes = pull(.data = control_val_data, var = formula.tools::get.vars(model)[1])
      
      treated_val_preds_trans_scale = sapply(X = 1:nrow(treated_train_fits),
                                             FUN = function(x) { 
                                               as.numeric(treated_val_mod_mats[["mod_mat"]][[x]] %*% coef(treated_train_fits[["mod"]][[x]])) })
      control_val_preds_trans_scale = sapply(X = 1:nrow(control_train_fits), 
                                             FUN = function(x) { 
                                               as.numeric(control_val_mod_mats[["mod_mat"]][[x]] %*% coef(control_train_fits[["mod"]][[x]])) })
      
      if(diff_k == 0){ 
        ## Use exp() to transform predictions back to original scale
        treated_val_preds = exp(treated_val_preds_trans_scale)
        control_val_preds = exp(control_val_preds_trans_scale)
        
      } else{
        ## Use exp() to transform predictions back to original scale
        treated_val_preds = exp((treated_val_preds_trans_scale) +
                                  pull(.data = filter(.data = treated_train_data, {{time_var}} == (val_time - diff_k)),
                                       var = formula.tools::get.vars(model)[1]))
        control_val_preds = exp((control_val_preds_trans_scale) +
                                  pull(.data = filter(.data = control_train_data, {{time_var}} == (val_time - diff_k)),
                                       var = formula.tools::get.vars(model)[1]))
      }
      
    }
  
  
  
  return(list("avg_treat_outcome" = mean(treated_val_outcomes),
              "avg_treat_pred" = mean(treated_val_preds),
              "avg_treat_pred_error" = mean(treated_val_outcomes) - mean(treated_val_preds),
              "avg_control_outcome" = mean(control_val_outcomes),
              "avg_control_pred" = mean(control_val_preds),
              "avg_control_pred_error" = mean(control_val_outcomes) - mean(control_val_preds),
              "avg_diff_pred_errors" = mean(treated_val_outcomes) - mean(treated_val_preds) - 
                (mean(control_val_outcomes) - mean(control_val_preds))))
  
}