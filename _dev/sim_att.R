sim_att <- function(fits, post_period, nsim = 200, cl = NULL, verbose = TRUE) {
    #In each simulation: compute prediction errors for each validation period using each model, compute differential average prediction errors
    
    #Compute original prediction errors
    time_var <- attr(fits, "time_var")
    data <- fits$data
    
    chk::chk_not_missing(post_period, "`post_period`")
    chk::chk_number(post_period)
    chk::chk_subset(post_period, data[[time_var]])
    
    unit_var <- attr(fits, "unit_var")
    group_var <- attr(fits, "group_var")
    group_levels <- sort(unique(data[[group_var]]))
    
    out_mat <- out_mat0 <- matrix(
        NA,
        nrow = length(fits$val_times),
        ncol = length(fits$models),
        dimnames = list(fits$val_times,
                        paste0("model_", seq_along(fits$models)))
    )
    
    observed_means <- predicted_means <- vector("list", nrow(fits$grid))
    
    for (i in seq_along(observed_means)) {
        val_time_ind <- fits$grid$val_time[i]
        mod <- fits$grid$model[i]
        fit <- fits$fits[[i]]
        
        val_data <- data[data[[time_var]] == fits$val_times[val_time_ind],, drop = FALSE]
        
        y <- model.response(model.frame(fit$formula, data = val_data))
        p <- predict(fit, newdata = val_data, type = "response")
        
        if (fits$models[[mod]]$log) {
            y <- exp(y)
            p <- exp(p)
        }
        
        observed_means[[i]] <- setNames(
            vapply(group_levels, function(g) mean(y[val_data[[group_var]] == g]), numeric(1L)),
            group_levels
        )
        
        predicted_means[[i]] <- setNames(
            vapply(group_levels, function(g) mean(p[val_data[[group_var]] == g]), numeric(1L)),
            group_levels
        )
        
        pred_error <- (observed_means[[i]]["1"] - predicted_means[[i]]["1"]) -
            (observed_means[[i]]["0"] - predicted_means[[i]]["0"])
        
        out_mat0[val_time_ind, mod] <- pred_error
    }
    
    #Select best based prediction errors
    worst_pred_within_model <- apply(abs(out_mat0), 2, max)
    best_model <- which.min(worst_pred_within_model)
    best_worst_pred <- worst_pred_within_model[best_model]
    
    #Compute ATT from best model
    model <- fits$models[[best_model]]
    post_fit <- fit_one_model(model$formula, data = data, family = model$family,
                         group_var = group_var, unit_var = unit_var,
                         time_var = time_var, val_time = post_period,
                         diff_k = model$diff_k, log = model$log)
    
    post_data <- data[data[[time_var]] == post_period,, drop = FALSE]
    
    y <- model.response(model.frame(post_fit$formula, data = post_data))
    p <- predict(post_fit, newdata = post_data, type = "response")
    
    if (model$log) {
        y <- exp(y)
        p <- exp(p)
    }
    
    post_observed_means <- setNames(
        vapply(group_levels, function(g) mean(y[post_data[[group_var]] == g]), numeric(1L)),
        group_levels
    )
    
    post_predicted_means <- setNames(
        vapply(group_levels, function(g) mean(p[post_data[[group_var]] == g]), numeric(1L)),
        group_levels
    )
    
    att <- (post_observed_means["1"] - post_predicted_means["1"]) -
        (post_observed_means["0"] - post_predicted_means["0"])
    
    
    #Draw parameters
    
    sim_coefs <- MASS::mvrnorm(nsim,
                               mu = unlist(fits$coefs),
                               Sigma = fits$vcov)
    
    # Compute prediction errors for each model for each validation period for each simulation
    
    chk::chk_flag(verbose)
    opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
    on.exit(pbapply::pboptions(opb))
    
    sim_list <- pbapply::pblapply(seq_len(nsim), function(s) {
        
        coefs <- sim_coefs[s,]
        
        for (i in seq_len(nrow(fits$grid))) {
            val_time_ind <- fits$grid$val_time[i]
            mod <- fits$grid$model[i]
            fit <- fits$fits[[i]]
            coefs_ind <- sum(lengths(fits$coefs[seq_along(fits$coefs) < i])) + seq_along(fits$coefs[[i]])
            
            #Compute pred error
            ##Create validation data
            val_data <- data[data[[time_var]] == fits$val_times[val_time_ind],, drop = FALSE]
            
            ##Set simulated coefficient
            fit <- marginaleffects::set_coef(fit, coefs[coefs_ind])
            
            ##Generate predictions on validation data
            p <- predict(fit, newdata = val_data, type = "response")
            
            #Unlog if outcome is logged to keep on original scale
            if (fits$models[[mod]]$log) {
                p <- exp(p)
            }
            
            predicted_means_s_i <- setNames(
                vapply(group_levels, function(g) mean(p[val_data[[group_var]] == g]), numeric(1L)),
                group_levels
            )
            
            pred_error <- (observed_means[[i]]["1"] - predicted_means_s_i["1"]) -
                (observed_means[[i]]["0"] - predicted_means_s_i["0"])
            
            out_mat[val_time_ind, mod] <- pred_error
        }
        
        #Select best based prediction errors
        worst_pred_within_model <- apply(abs(out_mat), 2, max)
        best_model <- which.min(worst_pred_within_model)
        best_worst_pred <- worst_pred_within_model[best_model]
        
        #Compute ATT from best model
        model <- fits$models[[best_model]]
        post_fit <- fit_one_model(model$formula, data = data, family = model$family,
                                  group_var = group_var, unit_var = unit_var,
                                  time_var = time_var, val_time = post_period,
                                  diff_k = model$diff_k, log = model$log)
        
        post_data <- data[data[[time_var]] == post_period,, drop = FALSE]
        
        y <- model.response(model.frame(post_fit$formula, data = post_data))
        p <- predict(post_fit, newdata = post_data, type = "response")
        
        if (model$log) {
            y <- exp(y)
            p <- exp(p)
        }
        
        post_observed_means <- setNames(
            vapply(group_levels, function(g) mean(y[post_data[[group_var]] == g]), numeric(1L)),
            group_levels
        )
        
        post_predicted_means <- setNames(
            vapply(group_levels, function(g) mean(p[post_data[[group_var]] == g]), numeric(1L)),
            group_levels
        )
        
        att <- (post_observed_means["1"] - post_predicted_means["1"]) -
            (post_observed_means["0"] - post_predicted_means["0"])
        
    }, cl = cl)
    
    sim_out <- simplify2array(sim_list, except = NULL)
    
    #Select best from simulated values
    
    
    
    out <- list(
        att = att,
        best_model = best_model,
        prediction_errors = out_mat0,
        sim_att = sim_att, #vector of att values from simulation
        sim_best_models = sim_best_models, #vector of index of best model in each simulation,
        sim_prediction_erros = sim_out #array of prediction errors in each simulation
    )
    
    out
}