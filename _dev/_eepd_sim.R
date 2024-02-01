.eepd_sim <- function(formula, data, family = gaussian(), vcov = NULL, 
                     group_var, unit_var, time_var, val_time, diff_k = 0,
                     log = FALSE, n_sim = 500, verbose = n_sim > 50) {
    
    # Argument checks
    chk::chk_not_missing(formula, "`formula`")
    chk::chk_is(formula, "formula")
    
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
    
    chk::chk_not_missing(val_time, "`val_time`")
    chk::chk_number(val_time)
    chk::chk_subset(val_time, data[[time_var]])
    
    chk::chk_count(diff_k)
    #need check to make sure diff_k isn't too high
    
    chk::chk_flag(log)
    
    chk::chk_count(n_sim)
    chk::chk_gt(n_sim, 0)
    
    chk::chk_flag(verbose)
    
    data[[unit_var]] <- factor(data[[unit_var]])
    
    data <- data[order(data[[unit_var]], data[[time_var]]),, drop = FALSE]
    
    # Create log and lagged variables
    outcome_name <- as.character(formula[[2]])
    
    formula_original <- formula
    
    # Log outcome if requested
    if (log) {
        outcome <- model.response(model.frame(formula, data = data))
        if (min(outcome) <= 0) {
            chk::err("`log` cannot be `TRUE` when the outcome takes on values of 0 or lower")
        }
        formula <- update(formula, log(.) ~ .)
    }
    
    # Add offset for lagged outcome if requested
    if (diff_k > 0) {
        lag_outcome_name <- sprintf("%s_lag_%s", outcome_name, diff_k)
        if (lag_outcome_name %in% all.vars(formula)) {
            chk::wrn(sprintf("the variable named %s will be replaced. Give this variable a different name before running",
                             diff_outcome_name))
        }
        
        for (i in seq_len(diff_k)) {
            lag_i <- data[[outcome_name]]
            is.na(lag_i)[] <- TRUE
            
            for (u in levels(data[[unit_var]])) {
                beg <- seq_len(i)
                end <- sum(data[[unit_var]] == u) + 1 - seq_len(i)
                lag_i[data[[unit_var]] == u][-beg] <- data[[outcome_name]][data[[unit_var]] == u][-end]
            }
            
            data[[lag_outcome_name]] <- {
                if (log) log(lag_i)
                else lag_i
            }
            
            formula <- update(formula, sprintf(". ~ . + offset(%s)", lag_outcome_name))
        }
    }
    
    train_data <- data[data[[time_var]] < val_time,, drop = FALSE]
    
    val_data <- data[data[[time_var]] == val_time,, drop = FALSE]
    
    train_data_1 <- droplevels(train_data[train_data[[group_var]] == 1,, drop = FALSE])
    train_data_0 <- droplevels(train_data[train_data[[group_var]] == 0,, drop = FALSE])
    
    val_data_1 <- droplevels(val_data[val_data[[group_var]] == 1,, drop = FALSE])
    val_data_0 <- droplevels(val_data[val_data[[group_var]] == 0,, drop = FALSE])
    
    val_data_1_by_unit <- split(val_data_1, val_data_1[[unit_var]])
    val_data_0_by_unit <- split(val_data_0, val_data_0[[unit_var]])
    
    if (is.character(family) && length(family) == 1 && family %in% c("negbin", "negative.binomial")) {
        fit_fun <- function(data, formula, family) {
            MASS::glm.nb(formula, data = data)
        }
    }
    else {
        fit_fun <- function(data, formula, family) {
            stats::glm(formula, data = data, family = family)
        }
    }
    
    #Model fits; one for each unit
    train_fits_1 <- lapply(split(train_data_1, train_data_1[[unit_var]]), fit_fun,
                           formula = formula, family = family)
    
    train_fits_0 <- lapply(split(train_data_0, train_data_0[[unit_var]]), fit_fun,
                           formula = formula, family = family)
    
    val_outcomes_1 <- val_data_1[[outcome_name]]
    val_outcomes_0 <- val_data_0[[outcome_name]]
    
    # Predictions for validation set for each unit
    val_preds_1 <- vapply(seq_along(train_fits_1), function(i) {
        marginaleffects::get_predict(train_fits_1[[i]],
                                     newdata = val_data_1_by_unit[[i]])[["estimate"]]
    }, numeric(1L))
    
    val_preds_0 <- vapply(seq_along(train_fits_0), function(i) {
        marginaleffects::get_predict(train_fits_0[[i]],
                                     newdata = val_data_0_by_unit[[i]])[["estimate"]]
    }, numeric(1L))
    
    if (log) {
        # Return predictions to original scale
        val_preds_1 <- exp(val_preds_1)
        val_preds_0 <- exp(val_preds_0)
    }
    
    if (n_sim > 0) {
        # Sample coefficients, get predictions
        
        ## Set up progress bar
        opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
        on.exit(pbapply::pboptions(opb))
        
        pb <- pbapply::startpb(0, n_sim * (length(train_fits_1) + length(train_fits_0)))
        
        val_preds_1_sim <- do.call("rbind", lapply(seq_along(train_fits_1), function(i) {
            #Fit for unit i
            fit_i <- train_fits_1[[i]]
            
            b <- marginaleffects::get_coef(fit_i)
            S <- marginaleffects::get_vcov(fit_i, vcov = vcov)
            
            coef_sim <- as.matrix(MASS::mvrnorm(n_sim, mu = b, Sigma = S))
            
            sapply(1:n_sim, function(s) {
                #Predictions for unit i for simulation s
                fit_i_s <- marginaleffects::set_coef(fit_i, coef_sim[s,])
                p <- marginaleffects::get_predict(fit_i_s, newdata = val_data_1_by_unit[[i]])[["estimate"]]
                pbapply::setpb(pb, n_sim * (i - 1) + s)
                p
            })
        }))
        
        val_preds_0_sim <- do.call("rbind", lapply(seq_along(train_fits_0), function(i) {
            #Fit for unit i
            fit_i <- train_fits_0[[i]]
            
            b <- marginaleffects::get_coef(fit_i)
            S <- marginaleffects::get_vcov(fit_i, vcov = vcov)
            
            coef_sim <- as.matrix(MASS::mvrnorm(n_sim, mu = b, Sigma = S))
            
            sapply(1:n_sim, function(s) {
                #Predictions for unit i for simulation s
                fit_i_s <- marginaleffects::set_coef(fit_i, coef_sim[s,])
                p <- marginaleffects::get_predict(fit_i_s, newdata = val_data_0_by_unit[[i]])[["estimate"]]
                pbapply::setpb(pb, n_sim * length(train_fits_1) + n_sim * (i - 1) + s)
                p
            })
        }))
        
        if (log) {
            # Return predictions to original scale
            val_preds_1_sim <- exp(val_preds_1_sim)
            val_preds_0_sim <- exp(val_preds_0_sim)
        }
        
        pbapply::setpb(pb, n_sim * (length(train_fits_1) + length(train_fits_0)))
        pbapply::closepb(pb)
    }
    
    out <- list(
        est = list(
            avg_treat_outcome = unname(mean(val_outcomes_1)),
            avg_treat_pred = unname(mean(val_preds_1)),
            avg_treat_pred_error = unname(mean(val_outcomes_1) - mean(val_preds_1)),
            avg_control_outcome = unname(mean(val_outcomes_0)),
            avg_control_pred = unname(mean(val_preds_0)),
            avg_control_pred_error = unname(mean(val_outcomes_0) - mean(val_preds_0)),
            avg_diff_pred_errors = unname(mean(val_outcomes_1) - mean(val_preds_1) - 
                                              (mean(val_outcomes_0) - mean(val_preds_0)))
        ),
        sim = if (n_sim > 0) list(
            avg_treat_outcome = unname(rep(mean(val_outcomes_1), n_sim)),
            avg_treat_pred = unname(colMeans(val_preds_1_sim)),
            avg_treat_pred_error = unname(mean(val_outcomes_1) - colMeans(val_preds_1_sim)),
            avg_control_outcome = unname(rep(mean(val_outcomes_0), n_sim)),
            avg_control_pred = unname(colMeans(val_preds_0_sim)),
            avg_control_pred_error = unname(mean(val_outcomes_0) - colMeans(val_preds_0_sim)),
            avg_diff_pred_errors = unname(mean(val_outcomes_1) - colMeans(val_preds_1_sim) - 
                                              (mean(val_outcomes_0) - colMeans(val_preds_0_sim)))
        )
    )
    
    attr(out, "info") <- list(
        formula = formula_original,
        family = family,
        group_var = group_var,
        unit_var = unit_var,
        time_var = time_var,
        val_time = val_time,
        diff_k = diff_k,
        log = log
    )
    
    class(out) <- "eepd_sim"
    
    out
}