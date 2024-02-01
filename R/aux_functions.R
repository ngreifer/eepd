#Fast (weighted) mean, optionally with subset
.wtd_mean <- function(x, w = NULL, subset = NULL) {
    if (is.null(subset)) {
        if (is.null(w)) {
            sum(x) / length(x)
        }
        else {
            sum(x * w) / sum(w)
        }
    }
    else {
        x <- x[subset]
        if (is.null(w)) {
            sum(x) / length(x)
        }
        else {
            w <- w[subset]
            sum(x * w) / sum(w)
        }
    }
}

.wtd_sd <- function(x, w = NULL, subset = NULL) {
    if (is.null(subset)) {
        if (is.null(w)) {
            sqrt(sum((x - .wtd_mean(x))^2)/(length(x) - 1))
        }
        else {
            sum_w <- sum(w)
            sqrt((sum_w / (sum_w^2 - sum(w^2))) * sum(w * (x - .wtd_mean(x, w))^2))
        }
    }
    else {
        x <- x[subset]
        if (is.null(w)) {
            sqrt(sum((x - .wtd_mean(x))^2)/(length(x) - 1))
        }
        else {
            w <- w[subset]
            sum_w <- sum(w)
            sqrt((sum_w / (sum_w^2 - sum(w^2))) * sum(w * (x - .wtd_mean(x, w))^2))
        }
    }
}

#Binds together multiple smaller square matrices (e.g., vcovs) into a larger one
#with 0s in the empty spaces. Used to create large covariance matrix
.block_diagonal <- function(matlist) {
    dim1 <- sum(unlist(lapply(matlist, ncol)))
    
    out <- matrix(0, nrow = dim1, ncol = dim1)
    k <- 0
    for (v in matlist) {
        ind <- k + seq_len(ncol(v))
        out[ind, ind] <- v
        k <- k + ncol(v)
    }
    
    out
}

#Checks if a given family specification is okay
.okay_family <- function(family) {
    if (is.character(family)) {
        if (length(family) != 1 || anyNA(family)) return(FALSE)
        if (family %in% c("negbin", "negative.binomial", "Negative Binomial")) return(TRUE)
        family <- get(family, mode = "function", envir = parent.frame(2))
    }
    if (is.function(family)) {
        family <- family()
    }
    
    !is.null(family$family) && is.function(family$variance) &&
        is.function(family$linkinv)
}

#Joint covariance matrix of coefficients across multiple models. Requires same units in all models,
#use nonzero weights to subset (e.g., weights of 1 for present and 1e-8 for absent). Should give
#same results as M-estimation (HC0 vcov).
vcovSUEST <- function(fits) {
    coef_lengths <- lengths(lapply(fits, function(f) na.omit(coef(f))))
    l <- c(0, cumsum(coef_lengths))
    coef_inds <- lapply(seq_along(fits), function(i) seq_len(coef_lengths[i]) + l[i])
    
    inf_func <- lapply(fits, function(f) { 
        b <- .bread(f)
        ef <- sandwich::estfun(f)
        inf <- tcrossprod(b, ef)/nobs(f)
        inf[is.na(inf)] <- 0
        inf
    })
    
    #VCOV matrix to be returned
    V <- matrix(NA_real_, nrow = sum(coef_lengths), ncol = sum(coef_lengths))
    
    for (i in seq_along(fits)) {
        ind_i <- coef_inds[[i]]
        
        #Usual within-model HC0 vcov
        V[ind_i, ind_i] <- tcrossprod(inf_func[[i]], inf_func[[i]])
        
        for (j in seq_along(fits)[-seq_len(i)]) {
            ind_j <- coef_inds[[j]]
            
            #between-model vcov components
            V[ind_i, ind_j] <- tcrossprod(inf_func[[i]], inf_func[[j]])
            V[ind_j, ind_i] <- t(V[ind_i, ind_j])
        }
    }
    
    V
}

#Quickly get bread matrix; for non-lm and nonglm objects, uses sandwich::bread()
.bread <- function(x) {
    if (!class(x)[1] %in% c("lm", "glm")) {
        return(sandwich::bread(x))
    }
    
    p <- x$rank
    
    if (p == 0) {
        return(matrix(NA_real_, 0L, 0L))
    }
    
    Qr <- x$qr
    
    coef.p <- x$coefficients[Qr$pivot[1:p]]
    cov.unscaled <- chol2inv(Qr$qr[1:p, 1:p, drop = FALSE])
    dimnames(cov.unscaled) <- list(names(coef.p), names(coef.p))
    
    df <- p + x$df.residual
    
    out <- cov.unscaled * df
        
    if (class(x)[1] == "glm" && !substr(x$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) {
        ww <- weights(x, "working")
        wres <- as.vector(residuals(x, "working")) * ww
        dispersion <- sum(wres^2, na.rm = TRUE) / sum(ww, na.rm = TRUE)
        out <- out * dispersion
    }
    
    out
}



