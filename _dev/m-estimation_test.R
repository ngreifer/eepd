#M-estimation test
library(geex)
data("lalonde", package = "MatchIt")

# Two arbtirary fits
fit1 <- glm(re78 > 1000 ~ age + educ, data = lalonde, family = binomial)
fit2 <- lm(re78  ~ age + educ + re74, data = lalonde)

# General M-estimation for coefficients of a (weighted) GLM
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

# Apply M-estimation to dataset and models (takes a bit of time)
estimates <- m_estimate(
    estFUN = eefun,
    data = lalonde,
    root_control = setup_root_control(start = 1.1 * c(coef(fit1), coef(fit2))),
    outer_args = list(modlist = list(fit1, fit2)))

all.equal(geex::coef(estimates), c(coef(fit1), coef(fit2)))

#Compare estimates
vcov_m <- geex::vcov(estimates)          #From geex
vcov_s1 <- sandwich::vcovHC(fit1, "HC0") #From sandwich for fit1
vcov_s2 <- sandwich::vcovHC(fit2, "HC0") #From sandwich for fit2
vcov_su <- eepd:::vcovSUEST(list(fit1, fit2))   #From eepd for fit1

# M-estimaiton vs. eepd
all.equal(vcov_m, vcov_su, check.attributes = FALSE)

# sandwich vs. eepd
all.equal(vcov_s1, vcov_su[1:length(coef(fit1)),
                           1:length(coef(fit1))], check.attributes = FALSE)
all.equal(vcov_s2, vcov_su[-(1:length(coef(fit1))),
                           -(1:length(coef(fit1)))], check.attributes = FALSE)

#Note: to use different subsets for each estimation, give each unit a weight
#of 1, and set to 0 any unit not to be included in an estimates. See below.
#Note in the package we set weights to NA instead which avoids a spurious
#warning. Set weights to a small value (e.g., 1e-9) to get sandwich
#results to agree, too.

# Subset to non-white units
lalonde$w1 <- 1; lalonde$w1[lalonde$race == "white"] <- 0
fit1 <- glm(re78 > 1000 ~ age + educ, data = lalonde, family = binomial,
            weights = w1)

# Subset to married units
lalonde$w2 <- 1; lalonde$w2[lalonde$married == 0] <- 0
fit2 <- lm(re78  ~ age + educ + re74, data = lalonde,
           weights = w2)

# Apply M-estimation to dataset and models (takes a bit of time)
estimates <- m_estimate(
    estFUN = eefun,
    data = lalonde,
    root_control = setup_root_control(start = 1.1 * c(coef(fit1), coef(fit2))),
    outer_args = list(modlist = list(fit1, fit2)))

all.equal(geex::coef(estimates), c(coef(fit1), coef(fit2)))

#Compare estimates
vcov_m <- geex::vcov(estimates)          #From geex
vcov_su <- eepd:::vcovSUEST(list(fit1, fit2))   #From eepd for fit1

# M-estimaiton vs. eepd
all.equal(vcov_m, vcov_su, check.attributes = FALSE)

