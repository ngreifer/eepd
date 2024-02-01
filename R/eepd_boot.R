#' @title Bootstrap estimation of ATTs
#' 
#' @description `eepd_boot()` bootstraps the selection of optimal models and estimation of ATTs done by [eepd_sim()] in order to account for uncertainty in sampling from the population. Bootstrapping is done by [fwb::fwb()], which uses the fractional weighted bootstrap or the traditional bootstrap.
#' 
#' @inheritParams eepd_fit
#' @inheritParams fwb::fwb
#' @param models either an `eepd_models` object (the output of a call to [eepd_mod()]) or an `eepd_fits` object (the output of a call to [eepd_fit()]). If the latter, the arguments `data`, `group_var`, `unit_var`, `time_var`, `val_times`, and `post_time` should be left empty as they will be extracted from the supplied object.
#' @param nboot the number of bootstrap iterations to use; default is 999. More is better but takes longer.
#' @param boot_type string; the type of bootstrap to perform. See the `wtype` argument of [fwb::fwb()] for allowable options. The default is `"exp"`, which requests the fractional weighted bootstrap using weights drawn from an Exp(1) distribution. `"multinom"` requests the usual bootstrap, which can fail when key observation requird to fit certain models happen not to be selected into a given bootstrap sample.
#' @param nsim the number of simulation iteration to perform in each bootstrap sample. Default is 200. More is better but takes longer.
#' 
#' @returns
#' An `fwb` object containing the estimated ATTs in each bootstrap iteration. See [fwb::fwb()] for details. `summary()`, `plot()`, and `print()` methods are available; see [fwb::summary.fwb()] and [fwb::plot.fwb()] for details.
#' 
#' @seealso [eepd_fit()]; [eepd_sim()]; [fwb::fwb()]
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models <- eepd_mod(list(crude_rate ~ 1,
#'                         crude_rate ~ year),
#'                    log = c(FALSE, TRUE))
#' models
#' 
#' # Fit the models to data; unit_var must be supplied for
#' # fixed effects
#' cl <- parallel::detectCores()
#' boot_out <- eepd_boot(models, data = ptpdata,
#'                       nboot = 99, nsim = 100,
#'                       group_var = "group",
#'                       time_var = "year",
#'                       val_times = 1999:2003,
#'                       post_time = 2008,
#'                       unit_var = "state",
#'                       cl = cl, verbose = TRUE)
#' 
#' summary(boot_out, ci.type = "perc")

#' @export 
eepd_boot <- function(models, data, nboot = 999, boot_type = getOption("fwb_wtype", "exp"), nsim = 200,
                      group_var, unit_var, time_var,
                      val_times, post_time, cl = NULL, verbose = FALSE) {
    # Argument checks
    mcall <- match.call()
    chk::chk_not_missing(models, "`models`")
    
    if (inherits(models, "eepd_models")) {
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
    }
    else if (inherits(models, "eepd_fits")) {
        chk::chk_missing(data, "`data`")
        chk::chk_missing(group_var, "`group_var`")
        chk::chk_missing(time_var, "`time_var`")
        chk::chk_missing(unit_var, "`unit_var`")
        chk::chk_missing(val_times, "`val_times`")
        chk::chk_missing(post_time, "`post_time`")
        
        data <- models$data
        group_var <- attr(models, "group_var")
        time_var <- attr(models, "time_var")
        unit_var <- attr(models, "unit_var")
        val_times <- models$val_times
        post_time <- models$post_time
        models <- models$models
    }
    else {
        chk::err("`models` must be an `eepd_models` object (the output of a call to `eepd_mod()`) or an `eepd_fits` object (the output of a call to `eepd_fit()`")
    }
    
    chk::chk_flag(verbose)
    
    .boot_fun <- function(.data, .weights) {
        
        fits <- .fit_models_internal(models, .data, .weights, group_var, unit_var, time_var,
                                     val_times, post_time)

        if (all(.weights == 1)) {
            sim_est <- eepd_sim(fits)
        }
        else {
            sim_est <- eepd_sim(fits, nsim)
        }
        
        c(ATT = unname(mean(sim_est$atts)))
    }
    
    boot_out <- fwb::fwb(data, .boot_fun, R = nboot, cluster = data[[unit_var]],
                         wtype = boot_type, verbose = verbose, cl = cl)
    
    class(boot_out) <- c("eepd_boot", class(boot_out))
    
    boot_out
}