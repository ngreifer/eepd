#' @title Fit models to datasets
#' 
#' @description `eepd_fit()` fits models to the pre-treatment and post-treatment models to obtain the joint distribution of model coefficients, which used by [eepd_sim()] to compute the average differential prediction errors and ATTs.
#' 
#' @param models an `eepd_models` object; the output of a call to [eepd_mod()].
#' @param data a dataset containing all the variables named in the supplied models (i.e., the outcome and any predictors) as well as any variable named below.
#' @param weights an optional vector of weights (e.g., sampling weights) used to fit weighted regression models.
#' @param group_var string; the name of the treatment variable in `data` defining the "to be treated" and "not to be treated" groups. The corresponding variable should take on values of 0 and 1 only.
#' @param unit_var string, optional; the name of the unit ID variable in `data` when any models suppiled to `models` contain unit fixed effects.
#' @param time_var string; the name of the variable in `data` containing the time variable.
#' @param val_times a numeric vector corresponding to the pre-treatment times that will be used as validation times when select the model with the optimal average expected prediction error.
#' @param post_time a number corresponding to the first post-treatment time period.
#' 
#' @returns
#' An `eepd_fits` object, which is a list containing the models supplied to `models`, a grid of all fitted models, a list of all model fit objects, a list of all estimated coefficients, the joint covariance of the coefficients, the dataset supplied to `data`, and other components supplied to `eepd_fit()`.
#' 
#' @details
#' `eepd_fit()` create a grid of all models and all time points and fits all corresponding models. For each validation time supplied to `val_times` and for the post-treatment time supplied to `post_time`, each model is fit using all previous times. For example, for a validation time of 5, a model is fit with data only from periods 1-4.
#' 
#' [glm()] or [MASS::glm.nb()] are used to fit the given models. The joint covariance matrix of all the coefficients is computed using the SUEST method described in Mize et al. (2019, p164), which is also used by the STATA command `suest`. This is equivalent to the covariance matrix computed by stacking the score equations for the models and fitting them using M-estimation and yields the equivalent of the HC0 covariance matrix for all within-model covariances.
#' 
#' @seealso [glm()]; [MASS::glm.nb()]; [eepd_sim()] to compute the ATT and average absolute differential prediction errors; [eepd_boot()] to bootstrap the models.
#' 
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models <- eepd_mod(list(crude_rate ~ 1,
#'                          crude_rate ~ year),
#'                     family = list("gaussian", "quasipoisson"),
#'                     lag = 0:1, fixef = TRUE)
#' models
#' 
#' # Fit the models to data; unit_var must be supplied for
#' # fixed effects
#' fits <- eepd_fit(models, data = ptpdata,
#'                  group_var = "group",
#'                  time_var = "year",
#'                  val_times = 1999:2007,
#'                  post_time = 2008,
#'                  unit_var = "state")


#' @export 
eepd_fit <- function(models, data, weights = NULL, group_var, time_var,
                     val_times, post_time, unit_var) {
    
    # Argument checks
    chk::chk_not_missing(models, "`models`")
    chk::chk_is(models, "eepd_models")
    
    chk::chk_not_missing(data, "`data`")
    chk::chk_data(data)
    
    # Process and order dataset
    data <- as.data.frame(data)
    
    if (is.null(rownames(data))) {
        rownames(data) <- seq_len(nrow(data))
    }
    
    chk::chk_not_missing(group_var, "`group_var`")
    chk::chk_string(group_var)
    chk::chk_subset(group_var, names(data))
    if (length(unique(data[[group_var]])) != 2) {
        chk::err("the grouping variable must have exactly 2 unique values")
    }
    data[[group_var]] <- factor(data[[group_var]])
    levels(data[[group_var]]) <- c("0", "1")
    
    chk::chk_not_missing(time_var, "`time_var`")
    chk::chk_string(time_var)
    chk::chk_subset(time_var, names(data))
    
    if (any(vapply(models, function(m) m$diff_k > 0 || m$lag > 0 || m$fixef, logical(1L)))) {
        chk::chk_not_missing(unit_var, "`unit_var`")
        chk::chk_string(unit_var)
        chk::chk_subset(unit_var, names(data))
        data[[unit_var]] <- factor(data[[unit_var]])
        
        data <- data[order(data[[unit_var]], data[[time_var]]),, drop = FALSE]
    }

    chk::chk_not_missing(val_times, "`val_times`")
    chk::chk_numeric(val_times)
    chk::chk_subset(val_times, data[[time_var]])
    
    chk::chk_not_missing(post_time, "`post_time`")
    chk::chk_number(post_time)
    chk::chk_subset(post_time, data[[time_var]])
    chk::chk_gt(post_time, val_times)
    
    .fit_models_internal(models = models,
                         data = data,
                         weights = weights,
                         group_var = group_var,
                         unit_var = unit_var,
                         time_var = time_var,
                         val_times = val_times,
                         post_time = post_time)
}

.fit_models_internal <- function(models, data, weights = NULL, group_var, unit_var, time_var,
                                 val_times, post_time) {
    
    #Note: data must be ordered correctly (by unit_var and time_var)
    if (is.null(weights)) {
        weights <- rep(1, nrow(data))
    }
    
    times <- c(val_times, post_time)
    
    # Create grid of models to be fit
    grid <- expand.grid(time_ind = seq_along(times), model = seq_along(models))
    
    #Remove models with problematic lags
    first_time <- min(data[[time_var]])
    lags <- pmax(vapply(models, `[[`, numeric(1L), "lag"), vapply(models, `[[`, numeric(1L), "diff_k"))
    lag_too_much <- times[grid$time_ind] - lags <= first_time

    if (any(lag_too_much)) {
        # if (all(lag_too_much)) {
            chk::err("some models involve lags corresponding to a period prior to the earliest time. Decrease the `lag` or `diff_k` components of the supplied models or use later validation times")
        # }
        # chk::wrn("some models involve lags corresponding to a period prior to the earliest time and will be omitted")
        # grid <- grid[!lag_too_much,, drop = FALSE]
    }
    
    #Fit all estimates
    coefs <- fits <- vector("list", nrow(grid))
    
    for (i in seq_len(nrow(grid))) {
        model <- models[[grid$model[i]]]
        val_time <- times[grid$time_ind[i]]
        fit <- .fit_one_model(model$formula, data = data, weights = weights,
                              group_var = group_var, unit_var = unit_var,
                              time_var = time_var, val_time = val_time,
                              model = model)
        fit$call <- NULL
        fits[[i]] <- fit
        coefs[[i]] <- coef(fit)
        # vcovs[[i]] <- sandwich::vcovHC(fit, type = vcov) #assuming separate; need to add m-estimation option
    }
    
    #If bootstrapping, do above many times
    #If M-estimation, compute giant vcov from fits
    
    # vcov <- .block_diagonal(vcovs) #assuming separate; need to add m-estimation option
    vcov <- vcovSUEST(fits) #joint distribution of coefs
    
    out <- list(models = models,
                val_times = val_times,
                post_time = post_time,
                grid = grid,
                fits = fits,
                coefs = coefs,
                vcov = vcov,
                data = data,
                weights = weights)
    
    attr(out, "time_var") <- time_var
    attr(out, "unit_var") <- unit_var
    attr(out, "group_var") <- group_var

    class(out) <- "eepd_fits"
    
    out
}

.fit_one_model <- function(formula, data, weights = NULL,
                           group_var, unit_var, time_var, val_time,
                           model) {
    
    # Create log and lagged variables
    mod <- .modify_formula_and_data(model, data, group_var, unit_var)
    
    # Effectively subset without dropping any observations
    # tol <- 1e-9
    # weights[weights < tol | mod$data[[time_var]] >= val_time] <- tol
    is.na(weights)[mod$data[[time_var]] >= val_time] <- TRUE
    
    # Model fitting function; note: need do.call() to correctly process `weights`
    if (model$family$family == "Negative Binomial") {
        fit_fun <- function(formula, data, family, weights = NULL, subset = NULL) {
            do.call(MASS::glm.nb, list(formula, data = data,
                                       weights = weights,
                                       subset = subset,
                                       na.action = "na.exclude"))
        }
    }
    else {
        fit_fun <- function(formula, data, family, weights = NULL, subset = NULL) {
            do.call(stats::glm, list(formula, data = data,
                                     family = family,
                                     weights = weights,
                                     subset = subset,
                                     na.action = "na.exclude"))
        }
    }
    
    fit_fun(mod$formula,
            data = mod$data,
            family = model$family,
            weights = weights)
}

.modify_formula_and_data <- function(model, data, group_var, unit_var) {
    
    formula <- model$formula
    
    # Create log and lagged variables
    outcome_name <- as.character(formula[[2]])
    
    # Add interaction with group
    formula <- update(formula, sprintf(". ~ %s * (.)", group_var))
    
    # Log outcome if requested
    if (model$log) {
        outcome <- model.response(model.frame(formula, data = data))
        if (min(outcome) <= 0) {
            chk::err("`log` cannot be `TRUE` when the outcome takes on values of 0 or lower")
        }
        formula <- update(formula, log(.) ~ .)
    }
    
    # Add lagged outcome or offset thereof if requested
    if (model$lag > 0 || model$diff_k > 0) {
        # 
        # if (lag_outcome_name %in% all.vars(formula)) {
        #     chk::wrn(sprintf("the variable named %s will be replaced. Give this variable a different name before running",
        #                      lag_outcome_name))
        # }
        
        for (i in seq_len(max(model$diff_k, model$lag))) {
            if (i > model$lag && i != model$diff_k) next
            
            lag_outcome_name <- sprintf("%s_lag_%s", outcome_name, i)
            
            lag_i <- data[[outcome_name]]
            is.na(lag_i)[] <- TRUE
            
            beg <- seq_len(i)
            for (u in levels(data[[unit_var]])) {
                #We can lag here because data is ordered by time_var already
                #Note: assumes complete time series for each unit, uses previous value in dataset (ignoring actual
                #      value of time var)
                end <- sum(data[[unit_var]] == u) + 1 - seq_len(i)
                lag_i[data[[unit_var]] == u][-beg] <- data[[outcome_name]][data[[unit_var]] == u][-end]
            }
            
            data[[lag_outcome_name]] <- lag_i
            
            if (i <= model$lag) {
                #Add lag as predictor
                formula <- {
                    if (model$log) {
                        update(formula, sprintf(". ~ . + log(%s)", lag_outcome_name))
                    }
                    else {
                        update(formula, sprintf(". ~ . + %s", lag_outcome_name))
                    }
                }
            }
            else if (i == model$diff_k) {
                #Add lag as offset
                formula <- {
                    if (model$log) {
                        update(formula, sprintf(". ~ . + offset(log(%s))", lag_outcome_name))
                    }
                    else {
                        update(formula, sprintf(". ~ . + offset(%s)", lag_outcome_name))
                    }
                }
            }
        }
    }
    
    # Add unit fixed effects if requested, remove group var main effect
    if (model$fixef) {
        formula <- update(formula, sprintf(". ~ . + %s - %s", unit_var, group_var))
    }
    
    list(formula = formula, data = data)
}