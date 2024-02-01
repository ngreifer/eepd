eefun <- function(data, modlist) {
    Xlist <- lapply(modlist, model.matrix, data = data)
    frames <- lapply(modlist, model.frame, data = data)
    Ylist <- lapply(frames, model.response)
    Wlist <- lapply(frames, model.weights)
    
    indexlist <- lapply(seq_along(modlist), function(i) {
        if (i == 1) seq_len(ncol(Xlist[[i]]))
        else sum(lengths(Xlist[seq_len(i - 1)])) + seq_len(ncol(Xlist[[i]]))
    })
    
    function(theta){
        unlist(lapply(seq_along(modlist), function(i) {
            Yhat <- Xlist[[i]] %*% theta[indexlist[[i]]]
            
            w <- if (length(Wlist[[i]]) == 0) rep(1, length(Yhat)) else Wlist[[i]]
            
            # #Canonical link
            # apply(Xlist[[i]], 2, function(x) sum((Ylist[[i]] - family(modlist[[i]])$linkinv(Yhat)) * x))
            
            apply(Xlist[[i]], 2, function(x) {
                sum(w * family(modlist[[i]])$mu.eta(Yhat) * x * (Ylist[[i]] - family(modlist[[i]])$linkinv(Yhat)) /
                        family(modlist[[i]])$variance(family(modlist[[i]])$linkinv(Yhat)))
            })
            
        }))
    }
}

library(geex)
fit1 <- glm(re78 > 1000 ~ age + educ, data = lalonde, family = binomial, weights = ifelse(married == 0, 1, 1e-9))
fit2 <- glm(re78 > 1100 ~ age + educ, data = lalonde, family = binomial)
estimates <- m_estimate(
    estFUN = eefun,
    data = lalonde,
    compute_roots = FALSE, roots = c(coef(fit1), coef(fit2)),
    # root_control = setup_root_control(start = c(coef(fit1), coef(fit2))),
    deriv_control = setup_deriv_control(method = "simple", method.args = list(eps = 1e-8)),
    outer_args = list(modlist = list((fit1),
                                     (fit2))))

# Compare point estimates
coef(estimates) # from GEEX
vcov(estimates)

.vcovSUEST(list(fit1, fit2))
vcovSUEST(list(fit1, fit2))
vcovSUR::vcovSUR(list(fit1, fit2), "rowid", F)


vcovsur <- function(ests) {
    inf_list <- lapply(ests, function(est) {
        X <- stats::model.matrix(est)
        if (any(alias <- is.na(stats::coef(est)))) 
            X <- X[, !alias, drop = FALSE]
        if (!is.null(w <- stats::weights(est))) 
            X <- X * sqrt(w)
        scores <- X * stats::resid(est)
        solve(crossprod(X), t(scores))
    })
    
    inf_funcs <- matrix(0, nrow = sum(unlist(lapply(inf_list, nrow))),
                        ncol = sum(unlist(lapply(inf_list, ncol))))
    k <- c(0, 0)
    for (i in seq_along(inf_list)) {
        inf_funcs[k[1] + seq_len(nrow(inf_list[[i]])),
                  k[2] + seq_len(ncol(inf_list[[i]]))] <- inf_list[[i]]
        
        k <- k + c(nrow(inf_list[[i]]), ncol(inf_list[[i]]))
    }
    
    # inf_funcs <- vcovSUR:::block_diagonal(inf_list)
    col_names <- unlist(lapply(inf_list, rownames))
    ks <- unlist(lapply(inf_list, nrow))
    
    cl <- unlist(lapply(ests, function(x) 1:nobs(x)))
    
    group_col_idx <- split(seq_along(cl), as.factor(cl))
    
    inf_funcs <- lapply(group_col_idx, function(idx) {
        rowSums(inf_funcs[, idx, drop = FALSE])
    })
    inf_funcs <- do.call("cbind", inf_funcs)
    
    vcov <- tcrossprod(inf_funcs)
    rownames(vcov) <- colnames(vcov) <- col_names
    
    vcov
}

fit1 <- glm(re78 > 1000 ~ age + married, data = lalonde_mis, family = binomial)
fit2 <- glm(re78 > 1000 ~ age + educ, data = lalonde_mis, family = binomial)
