#' @title Generate models used to fit outcomes
#' 
#' @description `eepd_mod()` generates a list of models characterized by a basic model formulas and other options (e.g., lags, families, etc.) that are supplied to [eepd_fit()], [eepd_sim()], or [eepd_boot()]. These values are completely crossed to create a grid of model specifications, and multiple sets of model specifications can be combined using `c()` (see Examples).
#' 
#' @param formula_list a list of model formulas with the outcome on the left side and predictions (or just an intercept) on the right side.
#' @param family a list of family specifications; see [family()] for allowable options. These iwll eventually be passed to [glm()] when fitting the models in [eepd_fit()], [eepd_sim()], or [eepd_boot()]. `"negbin"` can also be supplied to request a negative binomial model fit using [MASS::glm.nb()]. Default is `"gaussian"` to specify a linear model.
#' @param lag a vector of integers indicating the desired outcome lags to be used as predictors. For example, a `lag` value of 3 means the outcome lagged once, twice, and three times will be included as predictors. Default is 0 for no lags.
#' @param diff_k a vector of integers indicating the desired outcome lag to be used a an offset For example, a `diff_k` value of 3 means the third lag of the outcome will be included as an offset, equivalent to using the outcome minus its corresponding lag as the outcome of the corresponding model. Default is 0 for no lags. Any models with a `diff_k` value less than a `lag` value will be removed automatically.
#' @param log a logical vector indicating whether the outcome should be logged. Default is `FALSE` to use the original outcome.
#' @param fixef a logical vector indicating whether unit fixed effects should be included as predictors. Default is `FALSE` to omit unit fixed effects.
#' @param identiy_only_log `logical`; whether to omit any models in which `log` is `TRUE` but the link in the `family` specification corresponds to something other than `"identity"`. Default is `TRUE`, and this should probably not be changed.
#' 
#' @returns
#' An `eepd_models` object, which is a list containing the full cross (less any omitted combinations) of the model features specified in the arguments. These have a `print()` method and can be combined using `c()`.
#' 
#' @seealso [formula], [family]
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models1 <- eepd_mod(list(crude_rate ~ 1,
#'                          crude_rate ~ year),
#'                     family = list("gaussian", "quasipoisson"),
#'                     lag = 0:1, fixef = TRUE)
#' models1
#' 
#' # Add a single other model with a square time trend
#' models2 <- eepd_mod(crude_rate ~ poly(year, 2),
#'                     family = "gaussian",
#'                     fixef = FALSE)
#' models2
#' 
#' (models <- c(models1, models2))
#' 
#' # Remove a model
#' models[[4]] <- NULL
#' models


#' @export 
eepd_mod <- function(formula_list, family = "gaussian", lag = 0, diff_k = 0,
                     log = FALSE, fixef = FALSE, identiy_only_log = TRUE) {
    # Check arguments
    
    ## Check formula_list
    chk::chk_not_missing(formula_list, "`formula_list`")
    
    if (inherits(formula_list, "formula")) {
        formula_list <- list(formula_list)
    }
    else if (!is.list(formula_list) || !all(vapply(formula_list, inherits, logical(1L), "formula"))) {
        chk::err("`formula_list` must be a list of model formulas")
    }
    
    formula_list <- unique(formula_list)
    
    if (any(lengths(formula_list) < 3)) {
        chk::err("all formulas in `formula_list` must have left-hand-side (outcome) variable")
    }
    
    ## Check family
    if (.okay_family(family)) {
        family <- list(family)
    }
    else if (!is.list(family) ||
             !all(vapply(family, .okay_family, logical(1L)))) {
        chk::err("`family` must be a list of model families")
    }
    
    family <- unique(family)
    
    family <- lapply(family, function(f) {
        if (is.character(f)) {
            if (f %in% c("negbin", "negative.binomial", "Negative Binomial")) {
                return(list(family = "Negative Binomial", link = "log"))
            }
            
            f <- get(f, mode = "function", envir = parent.frame(2))
        }
        
        if (is.function(f)) {
            f <- f()
        }
        
        f
    })
    
    ## Check lag
    chk::chk_not_any_na(lag)
    chk::chk_whole_numeric(lag)
    chk::chk_gte(lag, 0)
    
    lag <- unique(lag)
    
    ## Check diff_k
    chk::chk_not_any_na(diff_k)
    chk::chk_whole_numeric(diff_k)
    chk::chk_gte(diff_k, 0)
    
    diff_k <- unique(diff_k)
    
    if (any(diff_k > 0)) {
        if (max(diff_k) <= min(lag))
            chk::wrn("`diff_k` will be ignored because all supplied values are less than the smallest value supplid to `lag`")
    }
    
    # Check log
    chk::chk_not_any_na(log)
    chk::chk_logical(log)
    
    log <- unique(log)
    
    # Check log
    chk::chk_not_any_na(fixef)
    chk::chk_logical(fixef)
    
    fixef <- unique(fixef)
    
    grid <- expand.grid(formula = seq_along(formula_list),
                        family = seq_along(family),
                        lag = seq_along(lag),
                        diff_k = seq_along(diff_k),
                        log = seq_along(log),
                        fixef = seq_along(fixef))
    
    if (any(log)) {
        chk::chk_flag(identiy_only_log)
        
        ## Remove all combinations that involve log = TRUE with a non-identity link (probably invalid)
        if (identiy_only_log) {
            identity_links <- vapply(family, function(f) {
                f[["link"]] == "identity"
            }, logical(1L))
            
            if (any(!identity_links)) {
                grid <- grid[grid$family %in% which(identity_links) | !log[grid$log],, drop = FALSE]
            }
        }
    }
    
    if (any(diff_k > 0) && any(lag > 0)) {
        #Drop combinations where diff_k is less than lag because those are equivalent to having no diff_k
        grid <- grid[diff_k[grid$diff_k] == 0 | diff_k[grid$diff_k] > lag[grid$lag],, drop = FALSE]
    }
    
    out <- lapply(seq_len(nrow(grid)), function(i) {
        list(formula = formula_list[[grid$formula[i]]],
             family = family[[grid$family[i]]],
             lag = lag[[grid$lag[i]]],
             diff_k = diff_k[[grid$diff_k[i]]],
             log = log[[grid$log[i]]],
             fixef = fixef[[grid$fixef[i]]])
    })
    
    class(out) <- "eepd_models"
    
    out
}

#' @exportS3Method print eepd_models
print.eepd_models <- function(x, ...) {
    for (i in seq_along(x)) {
        cat(sprintf("- Model %s:\n", i))
        cat(deparse1(x[[i]]$formula), "\n", sep = "")
        cat(sprintf("family: %s(link = '%s')\n",
                    x[[i]]$family$family,
                    x[[i]]$family$link))
        cat(sprintf("outcome lag: %s\n", if (x[[i]]$lag == 0) "none" else paste(seq_len(x[[i]]$lag), collapse = ", ")))
        cat(sprintf("outcome diff: %s\n", if (x[[i]]$diff_k == 0) "none" else x[[i]]$diff_k))
        cat(sprintf("log outcome: %s\n", if (x[[i]]$log) "yes" else "no"))
        cat(sprintf("unit fixed effects: %s\n", if (x[[i]]$fixef) "yes" else "no"))
        if (i != length(x)) cat("\n")
    }
    
    invisible(x)
}

#' @exportS3Method c eepd_models
c.eepd_models <- function(..., recursive = TRUE) {
    out <- NextMethod("c")
    
    out <- unique(out)
    
    class(out) <- "eepd_models"
    
    out
}