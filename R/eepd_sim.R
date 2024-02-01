#' @title Estimate ATTS from models fits
#' @name eepd_sim
#' 
#' @description `eepd_sim()` computes the ATTs from the models previously fit by [eepd_fit()], choosing the optimal one by minimizing the largest absolute average prediction error across validation times. Optionally, this process can be simulated to arrive at a distribution of ATTs that accounts for the uncertainty in selecting the optimal model. `plot()` plotst he resulting ATT(s).
#' 
#' @param fits an `eepd_fits` object; the output of a call to [eepd_fit()].
#' @param nsim the number of simulation iterations to run. Set to 0 (the default) to use the original model coefficients instead of simulating. See Details.
#' @param cl a cluster specification to allow for parallel computing. Passed to the `cl` argument of [pbapply::pblapply()]. Default is to use sequential evaluation on a single core.
#' @param verbose `logical`; when `nsim` is greater than 0, whether to print a progress bar. Default is `FALSE` to omit a progress bar.
#' @param x  an `eepd_sim` object; the output of a call to `eepd_sim()`.
#' @param stack `logical`; whether to produced a stacked density plot of the ATTs or multiple overlapping density plots when `nsim` is greater than 0.
#' @param palette string; the color palette to use when `nsim` is greater than 0. Passed to [ggplot2::scale_fill_brewer()].
#' @param ncol the number of columns in which to display the average absolute prediction errors when `nsim = 0`.
#' @param \dots ignored.
#' 
#' @returns
#' An `eepd_sim` object, which contains the absolute prediction errors for models in all validation periods in each simulation, the simulated ATTs for all models in each simulation, the optimal model in each simulation, and the chosen ATT in each simulation. When `nsim = 0`, these correspond to the  values computed using the original model coefficients.
#' 
#' `plot()` returns a `ggplot` object when `nsim = 0` and a `patchwork` object otherwise.
#' 
#' The `print()` method either displays the ATT and the model used to compute it (when `nsim = 0`) or the average of the simulated ATTs and a frequency table of the chosen models.
#' 
#' @details
#' `eepd_sim()` calculates the absolute average prediction errors for each validation year and the ATT in the post-treatment year specified in the supplied `eepd_fits` object using each model fit to the data by `eepd_fit()`. These quantities can be computed directly from the model coefficients and the data. The optimal model is selected as that which minimizes the absolute average prediction error across validation years, and the reported ATT is the ATT computed using that model. `plot()` can be used to display the absolute average prediction errors across years for each model. The black bar represents the year with the largest absolute average prediction error for the corresponding model.
#' 
#' Because there is uncertainty in which model is best, simulation can be used to obtain a distribution of optimal models and their corresponding ATTs. New sets of parameters are drawn from a multivariate normal distribution, the assumed "posterior" distribution of the model coefficients, and in simulation, absolute average prediction errors are computed, the optimal model is chosen, and the ATT from that model is computed. For each simulation, this yields an optimal model and its ATT. The final reported ATT is the average of the simulated ATTs, which may be computed from different models. To request this simulation, `nsim` should be set to a large number (e.g., 1000). `plot()` can be used to display the distribution of the ATTs; the plot contains both the distribution of the ATTs by model and the distribution of optimal models. The vertical dotted line represents the ATT reported by `print()`, which is the average of the simulated ATTs from just the optimal models chosen in each simulation (i.e., `mean(out$atts)` where `out` is the output of `eepd_sim()` with `nsim > 0`).
#' 
#' @seealso [eepd_fit()]; [eepd_boot()] to bootstrap the simulation process to correctly account for uncertainty sampling as well as in choosing the optimal model.
#' 
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models <- eepd_mod(list(crude_rate ~ 1,
#'                         crude_rate ~ year),
#'                    family = list("gaussian", "quasipoisson"),
#'                    lag = 0:1, fixef = TRUE)
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
#' 
#' # Choose the optimal model and compute its ATT
#' est <- eepd_sim(fits)
#' 
#' est
#' 
#' plot(est)
#' 
#' # Simulate average prediction errors and ATTs
#' est_sim <- eepd_sim(fits, nsim = 100)
#' 
#' est_sim
#' 
#' plot(est_sim)

#' @export
eepd_sim <- function(fits, nsim = 0, cl = NULL, verbose = FALSE) {
    #In each simulation: compute prediction errors for each validation period using each model, compute differential average prediction errors
    #Compute original prediction errors
    time_var <- attr(fits, "time_var")
    data <- fits$data
    
    unit_var <- attr(fits, "unit_var")
    group_var <- attr(fits, "group_var")
    group_levels <- sort(unique(data[[group_var]]))
    
    times <- c(fits$val_times, fits$post_time)
    
    weights <- fits$weights
    
    val_data <- val_weights <- vector("list", nrow(fits$grid))
    
    for (i in seq_len(nrow(fits$grid))) {
        mod <- fits$grid$model[i]
        d <- .modify_formula_and_data(fits$models[[mod]], data, group_var, unit_var)$data
        
        subset_i <- which(d[[time_var]] == times[fits$grid$time_ind[i]])
        
        val_data[[i]] <- d[subset_i,, drop = FALSE]
        val_weights[[i]] <- weights[subset_i]
        
        fits$coefs[[i]] <- na.omit(fits$coefs[[i]])
    }
    
    observed_means <- lapply(seq_len(nrow(fits$grid)), function(i) {
        mod <- fits$grid$model[i]
        fit <- fits$fits[[i]]
        
        y <- model.response(model.frame(update(fit$formula, . ~ 1), data = val_data[[i]]))
        
        if (fits$models[[mod]]$log) {
            y <- exp(y)
        }
        
        setNames(
            vapply(group_levels, function(g) {
                .wtd_mean(y, val_weights[[i]], val_data[[i]][[group_var]] == g)
            }, numeric(1L)),
            group_levels
        )
    })
    
    if (nsim > 0) {
        #Draw parameters
        sim_coefs <- MASS::mvrnorm(nsim,
                                   mu = unlist(fits$coefs),
                                   Sigma = fits$vcov)
        do_sim <- TRUE
    }
    else {
        sim_coefs <- matrix(unlist(fits$coefs), nrow = 1)
        nsim <- 1
        do_sim <- FALSE
    }
    
    # Compute prediction errors for each model for each validation period for each simulation
    
    if (do_sim) {
        chk::chk_flag(verbose)
    }
    opb <- pbapply::pboptions(type = if (verbose && do_sim) "timer" else "none")
    on.exit(pbapply::pboptions(opb))
    
    coefs_inds <- lapply(seq_len(nrow(fits$grid)), function(i) {
        sum(lengths(fits$coefs[seq_along(fits$coefs) < i])) + seq_along(fits$coefs[[i]])
    })
    
    mat0 <- matrix(NA_real_, nrow = length(times), ncol = length(fits$models),
                   dimnames = list(times,
                                   paste0("model_", seq_along(fits$models))))
    
    #out_mat: all prediction errors; length(times) x length(models) x nsim
    out_mat <- simplify2array(pbapply::pblapply(seq_len(nsim), function(s) {
        
        mat <- mat0
        
        coefs <- sim_coefs[s,]
        
        for (i in seq_len(nrow(fits$grid))) {
            mod <- fits$grid$model[i]
            fit <- fits$fits[[i]]
            
            #Compute pred error
            
            ##Set simulated coefficient
            fit <- marginaleffects::set_coef(fit, coefs[coefs_inds[[i]]])
            
            ##Generate predictions on validation data
            p <- predict(fit, newdata = val_data[[i]], type = "response")
            
            #Unlog if outcome is logged to keep on original scale
            if (fits$models[[mod]]$log) {
                p <- exp(p)
            }
            
            predicted_means_s_i <- setNames(
                vapply(group_levels, function(g) {
                    .wtd_mean(p, val_weights[[i]], val_data[[i]][[group_var]] == g)
                }, numeric(1L)),
                group_levels
            )
            
            pred_error <- (observed_means[[i]]["1"] - predicted_means_s_i["1"]) -
                (observed_means[[i]]["0"] - predicted_means_s_i["0"])
            
            mat[fits$grid$time_ind[i], mod] <- pred_error
            
            # out_mat[fits$grid$time_ind[i], mod, s] <- pred_error
        }
        
        mat
    }, cl = cl))
    
    out_mat_att <- out_mat[times == fits$post_time,,, drop = FALSE]
    
    out_mat <- abs(out_mat[times != fits$post_time,,, drop = FALSE])
    
    optimal_models <- vapply(seq_len(nsim), function(s) {
        worst_pred_within_model <- {
            if (is.null(dim(out_mat[,, s]))) max(out_mat[,, s])
            else apply(out_mat[,, s], 2, max)
        }
        which.min(worst_pred_within_model)
    }, integer(1L))
    
    atts <- out_mat_att[cbind(1, optimal_models, seq_len(nsim))]
    
    out <- list(sim_pred_errors = out_mat,
                sim_atts = out_mat_att,
                optimal_models = optimal_models,
                atts = atts,
                models = fits$models)
    
    attr(out, "simulated") <- do_sim
    
    class(out) <- "eepd_sim"
    
    out
}

#' @exportS3Method print eepd_sim
print.eepd_sim <- function(x, digits = getOption("digits", 3), ...) {
    n_mods <- dim(x$sim_atts)[2]
    opt_mods <- factor(x$optimal_models,
                       levels = seq_len(n_mods),
                       labels = paste("Model", seq_len(n_mods)))
    
    if (attr(x, "simulated")) {
        cat("Simulation results:\n")
        cat(sprintf("- ATT: %s\n", format(mean(x$atts), digits = digits)))
        cat("- Optimal models:\n")
        print(table(opt_mods))
    }
    else {
        cat("Estimation results:\n")
        cat(sprintf("- ATT: %s\n", format(x$atts, digits = digits)))
        cat(sprintf("- Optimal model: %s\n", as.character(opt_mods)))
    }
    
    invisible(x)
}

#' @rdname eepd_sim
#' @exportS3Method plot eepd_sim
plot.eepd_sim <- function(x, stack = TRUE, palette = "YlGnBu", ncol = NULL, ...) {
    chk::chk_flag(stack)
    chk::chk_string(palette)
    
    atts <- x$atts
    n_mods <- dim(x$sim_atts)[2]
    opt_mods <- factor(x$optimal_models,
                       levels = seq_len(n_mods),
                       labels = paste("Model", seq_len(n_mods)))
    
    opt_mods <- droplevels(opt_mods)
    
    if (attr(x, "simulated")) {

        est <- mean(atts)
        
        # Augment ATTs with low-weighted units to account for small counts
        w <- rep(1, length(opt_mods))
        tol <- 8
        
        for (i in levels(opt_mods)) {
            num_i <- sum(opt_mods == i)
            if (num_i <= tol) {
                atts <- c(atts, rep(mean(atts[opt_mods == i]), tol - num_i))
                opt_mods <- c(opt_mods, factor(rep(i, tol - num_i)))
                w <- c(w, rep(1e-8, tol - num_i))
                w[opt_mods == i] <- w[opt_mods == i] * num_i / sum(w[opt_mods == i])
            }
        }
        
        # Select bandwidth using weighted ATTs
        bw <- 1.06 * .wtd_sd(atts, w) * sum(w)^(-1/5)
        
        # Density plot
        p1 <- ggplot() + 
            geom_density(aes(x = atts,
                             fill = opt_mods,
                             weight = w,
                             y = after_stat(count)/sum(after_stat(count))),
                         alpha = if (stack) .8 else .3,
                         position = if (stack) "stack" else "identity",
                         bw = bw) +
            geom_vline(xintercept = est, linetype = "dashed") +
            geom_hline(yintercept = 0) +
            labs(x = "ATT", y = "Density") +
            guides(fill = "none") +
            scale_fill_brewer(palette = palette) +
            theme_bw() +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid = element_blank())
        
        # Bar plot
        p2 <- ggplot() +
            geom_bar(aes(x = 0, fill = opt_mods, weight = w,
                         y = after_stat(count)/sum(after_stat(count)))) +
            coord_cartesian(expand = FALSE) +
            labs(y = "Proportion", fill = "Optimal Model") +
            scale_fill_brewer(palette = palette) +
            theme_bw() +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.border = element_rect(fill = NA, color = "black"))
        
        # Put them together
        patchwork::wrap_plots(p1, p2, guides = "collect",
                              widths = c(5, 1))
    }
    else {
        df <- expand.grid(time = seq_len(dim(x$sim_pred_errors)[1]),
                          model = seq_len(dim(x$sim_pred_errors)[2]))
        
        df$pred_error <- x$sim_pred_errors[cbind(as.matrix(df), 1)]
        df$time <- factor(dimnames(x$sim_pred_errors)[[1]][df$time])
        df$model <- factor(dimnames(x$sim_pred_errors)[[2]][df$model],
                           levels = dimnames(x$sim_pred_errors)[[2]],
                           labels = gsub("model_", "Model ", dimnames(x$sim_pred_errors)[[2]], fixed = TRUE))
        df$is_max <- factor("no", levels = c("no", "yes"))
        
        for (j in levels(df$model)) {
            df$is_max[df$model == j][which.max(df$pred_error[df$model == j])] <- "yes"
        }
        
        p <- ggplot(df) +
            geom_col(aes(x = .data$time, y = .data$pred_error, fill = .data$is_max),
                     width = .96) +
            geom_hline(yintercept = 0) +
            facet_wrap(vars(.data$model), ncol = ncol) +
            scale_fill_manual(values = c("yes" = "black", "no" = "gray70")) +
            labs(y = "Absolute Difference in Average Prediction Errors",
                 x = "Validation Year") +
            guides(fill = "none") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = .5))
        
        p
    }
}