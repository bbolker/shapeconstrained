library(scam)
library(glmmTMB)
library(marginaleffects)
library(ggplot2); theme_set(theme_bw())
library(purrr)
library(dplyr)
library(bbmle)
library(RTMB)
library(tmbstan)
rstan_options(auto_write=TRUE) ## threads_per_chain not active
options(mc.cores = 4)

## what for tidy predictions? ggpredict, marginaleffects, emmeans?
source("funs.R")

if (!interactive()) pdf("reedfrog.pdf")
## GOAL: fit Reed frog predation data, with CIs, using Holling type 2,
##  unrestricted GAM, shape-constrained GAM (both scam and RTMB), ...
##

## tests? parametric bootstrap, AIC, etc. ?
## (edf for RTMB fits?)

set.seed(101)

data("ReedfrogFuncresp", package = "emdbook")
dd <- ReedfrogFuncresp
ddx <- expand_bern(dd, response = "Killed", size = "Initial")

ddp0 <- data.frame(Initial = 1:100)
## works for scam/gam fits
pfun <- function(m) (c(marginaleffects::predictions(m, newdata = ddp0))
    |> as.data.frame()
    |> as_tibble()
    |> dplyr::select(prob = estimate, lwr = conf.low, upr = conf.high)
    |> mutate(Initial = ddp0$Initial, .before = 1)
)

## unconstrained GAM fit
m_gam_tp <- gam(cbind(Killed, Initial -Killed) ~ s(Initial, bs = "tp",
                                                          k = 8), data = dd,
                       family = binomial)
preds_gam_tp <- pfun(m_gam_tp)

## This doesn't work -- don't know if it's a thinko or just numerical nastiness ...
## holling  a*x/(1+a*h*x)
## initial slope = a; asymptote = 1/h ~ 0.5, 100
## holling prob = a/(1+a*h*x)
## holling via GLM; partial fractions!
## 1/prob = A + B/x = (Ax + B)/x ->
## prob = x/(Ax + B) = (1/A)/(1 + (B/A)*(1/x))
## a = 1/A; h = B

## 1/prob ~ b0 + b1/x
## prob = 1/(b0 + b1/x) = (1/b0)/(1/(1/b0) + (1/b0)*b1
if (FALSE) {
    ## don't want this lying around messing things up
    m_glm_holling <- glm(cbind(Killed, Initial -Killed) ~ I(1/Initial), data = dd,
                         family = binomial(link = "inverse"),
                         start = c(2, 0.02))
}
## coef(glm_holling)  ## bogus ...
## negative values for b1 ???
## plot(predict(glm_holling, newdata = data.frame(Initial = 1:100), type = "response"))


## fit Holling type 2 with mle2
m_mle2_holling <- bbmle::mle2(Killed ~ dbinom(prob = exp(loga)/(1+exp(loga)*exp(logh)*Initial),
                                                     size = Initial),
                                     start = list(loga = log(0.5), logh = log(0.01)),
                                     data = dd)

## predict number killed, divide by value
p0 <- predict(m_mle2_holling, newdata = list(Initial = 1:100))/(1:100)
rpars <- MASS::mvrnorm(1000, mu = coef(m_mle2_holling), Sigma = vcov(m_mle2_holling))
n_initial <- 1:100
preds <- apply(rpars, 1,
               function(x) { a <- exp(x[1]); h <- exp(x[2]); a/(1+a*h*n_initial)})
## matplot(preds, type = "l", col = adjustcolor("black", alpha = 0.3), lty = 1)
ci <- t(apply(preds, 1, quantile, c(0.025, 0.975)))


preds_mle2_holling <- data.frame(Initial = 1:100, prob = p0, lwr = ci[,1], upr = ci[,2])

## fit Holling type 2 with RTMB
m_RTMB_holling <- fit_RTMB_holling2(dd)
preds_RTMB_holling <- predict_RTMB_holling2(dd, data.frame(Initial = 1:100),
                                         m_RTMB_holling$fit$par) |>
    dplyr::select(Initial, prob, lwr, upr)

## fit SCAM with RTMB (Laplace approx doesn't work, need random = NULL)
m_RTMB_mpd <- fit_mpd_fun(data = dd, response = "Killed",
                          size = dd$Initial, xvar = "Initial",
                          family = "binomial", random = NULL)


m_RTMB2_mpd <- fit_mpd_fun(data = dd, response = "Killed",
                          size = dd$Initial, xvar = "Initial",
                          family = "binomial", random = "b1",
                          inner.control = list(smartsearch=FALSE, maxit =1),
                          opt  = "BFGS")

m_RTMB_mpd$fit
m_RTMB2_mpd$fit


## don't go outside original range (i.e initial<5), technical issues with outer.ok in splineDesign ...
dnew <- data.frame(Initial = 5:100)

rf_predfun <- function(fit, newdata = dnew, olddata = dd, ci_level = 0.95, random = NULL, extra_pars = NULL, ...) {
    n_new <- nrow(newdata)
    ddp <- data.frame(Initial = c(olddata$Initial, newdata$Initial),
                      Killed = c(olddata$Initial, rep(NA_integer_, n_new)))
    k <- data.frame(Initial = smoothCon(s(Initial, bs="mpd"), data = olddata, absorb.cons = TRUE)[[1]]$knots)
    ee <- fit$obj$env
    pp <- c(split(unname(ee$last.par.best), names(ee$last.par.best)), extra_pars)
    preds0 <- fit_mpd_fun(data = ddp, response = "Killed",
                          size = ddp$Initial,
                          xvar = "Initial",
                          family = "binomial",
                          random = random,
                          knots = k,
                          predict = TRUE,
                          parms = pp, ...)
    qq <- qnorm((1+ci_level)/2)
    ret <- (preds0
        |> filter(nm == "eta")
        |> slice_tail(n = n_new)
        |> transmute(Initial = newdata$Initial, prob = plogis(value), lwr = plogis(value-qq*sd),
                     upr = plogis(value+qq*sd))
    )
    return(ret)
}

preds_RTMB_mpd <- rf_predfun(m_RTMB_mpd)
## FIXME: still need to sort out warning messages
preds_RTMB2_mpd <- rf_predfun(m_RTMB2_mpd, random = "b1",  inner.control = list(smartsearch=FALSE, maxit =1))

m_scam_mpd <- scam(Killed ~ s(Initial, bs = "mpd"), data = ddx, family = binomial)

preds_scam_mpd <- pfun(m_scam_mpd)

(all_models <- ls(pattern="^m_"))
(all_preds <- ls(pattern="^preds_"))

## NOT 'preds_' (don't pollute namespace)
pred_frame <- (mget(all_preds)
    |> setNames(all_preds)
    |> bind_rows(.id = "model")
    |> mutate(across(model, ~gsub("preds_", "", .)))
)

pred_plot <- function(var, data = pred_frame) {
    var <- enquo(var)
    gg0 <- ggplot(data, aes(Initial, prob)) +
        geom_line(aes(colour = !!var)) +
        geom_ribbon(aes(ymin = lwr, ymax = upr, fill = !!var), colour = NA, alpha = 0.5) +
        expand_limits(y=0) +
        geom_point(data=dd, aes(y = Killed/Initial, size = Killed), alpha = 0.5) +
        facet_wrap(vars(!!var)) +
        theme(legend.position = "none")
    return(gg0)
}

print(pred_plot(model))

## compute log-likelihoods?
scam_df <- attr(logLik(m_scam_mpd), "df")
pp <- drop(predict(m_scam_mpd, newdata = dd, type = "response"))
scam_nll <- with(dd, -sum(dbinom(Killed, size = Initial, prob = pp, log = TRUE)))
scam_AIC <- 2*(scam_nll + scam_df)

gam_df <- attr(gll <- logLik(m_gam_tp), "df")
gam_nll <- -1*c(gll)
gam_AIC <- AIC(m_gam_tp)

mle2_nll <- -1*(mll <- logLik(m_mle2_holling))
mle2_df <- attr(mll, "df")
mle2_AIC <- AIC(m_mle2_holling)

with(dd, -sum(dbinom(Killed, size = Initial, prob = m_RTMB_mpd$mu, log = TRUE)))

rf_aictab <- tibble(
    model = c("scam/mpd", "gam/tp", "mle2/holling"),
    AIC = c(scam_AIC, gam_AIC, mle2_AIC),
    nll = c(scam_nll, gam_nll, mle2_nll),
    df = c(scam_df, gam_df, mle2_df)
) |> mutate(across(c(AIC, nll), ~ . - min(.)))

save("dd", "pred_plot", "pred_frame", "rf_aictab", file = "reedfrog_stuff.rda")
## scam CIs are not monotonic??

## investigate 'oversmoothing' of RTMB ...
## FIXME: clean up/integrate with rf_predfun above, OR get rid of it ?

if (FALSE) {
    ## w: working weights (precisions of params?)
    
     # calculating tr(A)...
     I.plus <- rep(1,nobs)   # define diagonal elements of the matrix I^{+}
     I.plus[w<0] <- -1
     L <- c(1/alpha)    # define diagonal elements of L=diag(1/alpha)
     ## NOTE PKt is O(np^2) and not needed --- can compute trA as side effect of gradiant
     KtILQ1R <- crossprod(L*I.plus*K,wX1) ## t(L*I.plus*K)%*%wX1 
     edf <- rowSums(P*t(KtILQ1R))
     trA <- sum(edf)
}

ddp <- data.frame(Initial = c(dd$Initial, 5:100),
                  Killed = c(dd$Initial, rep(NA_integer_, 96)))
k <- data.frame(Initial = smoothCon(s(Initial, bs="mpd"), data = dd, absorb.cons = TRUE)[[1]]$knots)
get_mpd_fix_preds <- function(log_smSD, ci_level = 0.95) {
    nsm <- sum(names(m_RTMB_mpd$fit$par) == "b1")
    fixparms <- list(
        b0 = 0,
        b1 = rep(0, nsm),
        ## log(unname(drop(m_scam_mpd$sp)))/2  ## drop extra dimensions etc
        log_smSD = log_smSD
    )
    m_RTMB_mpd_fix <- fit_mpd_fun(data = dd, response = "Killed",
                                  parms = fixparms,
                                  map = list(log_smSD = factor(NA_real_)),
                                  size = dd$Initial, xvar = "Initial",
                                  family = "binomial",
                                  random = NULL)
     preds1 <- fit_mpd_fun(data = ddp, response = "Killed",
                           size = ddp$Initial, xvar = "Initial",
                           family = "binomial", random = NULL,
                           knots = k, predict = TRUE,
                           parms = with(m_RTMB_mpd_fix$obj$env, parList(last.par.best)))
    qq <- qnorm((1+ci_level)/2)
    preds_RTMB_mpd_fix <- (preds1
        |> filter(nm == "eta")
        |> slice_tail(n = 96)
        |> transmute(Initial = 5:100, prob = plogis(value), lwr = plogis(value-qq*sd),
                     upr = plogis(value+qq*sd))
    )
    return(preds_RTMB_mpd_fix)
}


## pred breaks (no SDs) if log_smSD >= 2
sdvec <- c(1.5, 1, 0, -1, -2, -4)
pred_fix <- purrr::map_dfr(setNames(sdvec, sdvec),
               get_mpd_fix_preds,
               .id = "log_smSD") |> mutate(across(log_smSD, as.numeric))

print(pred_plot(log_smSD, pred_fix))

if (!interactive()) dev.off()

## not working (wrong number of knots)
if (FALSE) m_RTMB_mpd <- fit_mpd_fun(data = dd, response = "Killed",
                          size = dd$Initial, xvar = "Initial",
                          knots = data.frame(Initial = unique(dd$Initial)),
                          family = "binomial", random = NULL)


## debug(smooth.construct.mpd.smooth.spec)
## dk <- unique(dd["Initial"])
## nk <- nrow(dk)
## smoothCon(s(Initial, bs = "mpd", k=4), data = dd, absorb.cons = TRUE)


if (FALSE) {
    s1 <- tmbstan(m_RTMB_mpd$obj)
    library(shinystan)
    launch_shinystan(s1)
}

