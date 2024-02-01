.apply_models <- function(models, data, vcov = "HC2", group_var, unit_var, time_var,
                         val_times, n_sim = 500, cl = NULL, verbose = FALSE) {
    
    # Argument checks
    chk::chk_not_missing(models, "`models`")
    chk::chk_is(models, "eepd_models")
    
    chk::chk_not_missing(data, "`data`")
    chk::chk_data(data)
    
    chk::chk_not_missing(group_var, "`group_var`")
    chk::chk_string(group_var)
    chk::chk_subset(group_var, names(data))
    
    chk::chk_not_missing(unit_var, "`unit_var`")
    chk::chk_string(unit_var)
    chk::chk_subset(unit_var, names(data))
    
    chk::chk_not_missing(time_var, "`time_var`")
    chk::chk_string(time_var)
    chk::chk_subset(time_var, names(data))
    
    chk::chk_not_missing(val_times, "`val_times`")
    chk::chk_numeric(val_times)
    chk::chk_subset(val_times, data[[time_var]])
    
    chk::chk_count(diff_k)
    #need check to make sure diff_k isn't too high
    
    chk::chk_flag(log)
    
    chk::chk_count(n_sim)
    chk::chk_gt(n_sim, 0)
    
    chk::chk_flag(verbose)
    
    data[[unit_var]] <- factor(data[[unit_var]])
    
    data <- data[order(data[[unit_var]], data[[time_var]]),, drop = FALSE]
    
    grid <- expand.grid(val_time = seq_along(val_times), model = seq_along(models))
    
    fits <- pbapply::pblapply(seq_len(nrow(grid)), function(i) {
        val_time <- val_times[grid$val_time[[i]]]
        model <- models[[grid$model[[i]]]]
        eepd_sim(model$formula, data = data, family = model$family,
                 vcov = vcov, group_var = group_var, unit_var = unit_var,
                 time_var = time_var, val_time = val_time,
                 diff_k = model$diff_k, log = model$log, n_sim = n_sim, verbose = FALSE)
    })
    
    ests <- do.call("rbind", lapply(fits, function(f) as.data.frame(f$est)))
    
    sims <- do.call("rbind", lapply(fits, function(f) as.data.frame(f$est)))
    
}
