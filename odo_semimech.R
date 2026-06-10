mk_semimech_fun <- function(data, parms, random = "b1", silent = TRUE,
                       family = "gaussian", ...) {
    ## can't use %~% format if we want to add a penalty
    f <- function(parms) {
        getAll(data, parms)
        ## smoothing params positive for monotonicity
        b_ar <- exp(b_ar)
        b_ht <- exp(b_ht)
        ## exponentiate to maintain positivity
        ar <- exp(b0_ar + X_ar %*% b_ar)
        ht <- exp(b0_ht + X_ht %*% b_ht)
        mu <- 1/(1/ar + ht*initial) ## == ar/(1 + ar*ht*initial)
        eta <- qlogis(mu)
        nll <- 0
        for (i in 1:length(y)) {
            if (!is.na(y[i])) {
                nll <- nll -1*dbinom(y[i], prob = mu[i], size = initial[i], log = TRUE)
            }
        }
        pen <- (exp(-2*log_sm_ar) * (t(b_ar) %*% S_ar %*% b_ar) + 2*log_sm_ar)/2 +
               (exp(-2*log_sm_ht) * (t(b_ht) %*% S_ht %*% b_ht) + 2*log_sm_ht)/2
        REPORT(mu)
        ADREPORT(mu)
        REPORT(eta)
        ADREPORT(eta)
        return(nll + pen)
    }
    MakeADFun(f, parms, random=random, silent = silent, ...)
}

fit_semimech_fun <- function(data,
                        response = "y",
                        size = numeric(0),
                        parms = NULL,
                        knots = NULL,
                        predict = FALSE,
                        random = c("b_ar","b_ht"),
                        silent = TRUE,
                        opt = "nlminb",
                        se.fit = TRUE,
                        ...) {
    ## if predicting, make sure to pass old knots so basis is constructed properly
    sm_ar <- smoothCon(s(Size, bs = "cv"), data = data, absorb.cons = TRUE, knots = knots[[1]])[[1]]
    sm_ht <- smoothCon(s(Size, bs = "mpd"), data = data, absorb.cons = TRUE, knots = knots[[2]])[[1]]
    parms <- list(b0_ar = 0, b0_ht = 0,
                  b_ar = rep(0, ncol(sm_ar$X)),
                  b_ht = rep(0, ncol(sm_ht$X)),
                  log_sm_ar = 0, log_sm_ht = 0)
    data$y <- data$killed
    tmbdat <- c(as.list(data),
                list(size = size,
                     S_ar = sm_ar$S[[1]], X_ar = sm_ar$X,
                     S_ht = sm_ht$S[[1]], X_ht = sm_ht$X)
                )
    obj <- mk_semimech_fun(data = tmbdat, parms = parms,
                      random = random, 
                      silent = silent, ...)
    ## p0 <- parms[names(parms) != "b1"]
    ## optim(par = p0, fn = obj$fn, control = list(maxit = 2000))
    if (predict) {
        ## shouldn't need to map() b since we are using best-fit  if random = NULL ?
        ## if (!is.null(random)) parms <- parms[setdiff(names(parms), random)]
        obj$fn(unlist(parms))
        if (se.fit) {
            sdr <- sdreport(obj)
            return(with(sdr,
                        data.frame(nm = names(value), value, sd)))
        } else {
        }
    }
    res <- with(obj,
                switch(opt,
                       nlminb =  try(nlminb(par, fn, gr, control = list(eval.max = 1000, iter.max = 1000))),
                       BFGS = try(optim(par, fn, gr, method = "BFGS", control = list(maxit =1000))),
                       stop("unknown optimizer ", opt))
                )
    ret <- c(list(fit = res, obj = obj), obj$report())
    class(ret) <- c("myRTMB", "list")
    return(ret)
}
    
###
library(tidyverse); theme_set(theme_bw())
library(mgcv)
library(scam)
library(RTMB)

source("funs.R")

datfn <- "McCoy_response_surfaces_Gamboa.csv"

## odonates only
x <- (read.csv(datfn)
    |> transform(block = factor(block))
    |> subset(cohort == "single" & predtype == "odo",
              select = -c(cohort, predtype))
    |> droplevels()
    |> transform(csize = size - mean(size), prop = killed/initial)
    |> rename(Size = "size")
)

m_RTMB_sm <- fit_semimech_fun(data = x)
## m_RTMB_sm_optim <- fit_semimech_fun(data = x, opt = "BFGS")
## m_RTMB_sm_optim
## these fail
## fit_semimech_fun(data = x, inner.control = list(smartsearch=FALSE, maxit =1))
## fit_semimech_fun(data = x, inner.control = list(smartsearch=FALSE, maxit =1), opt = "BFGS")

## adapt rf_predfun guts
newdata <- expand.grid(Size = 8:50, initial = 0:100)
olddata <- x
init_dens <- "initial"
response <- "killed"
fit <- m_RTMB_sm
n_new <- nrow(newdata)

ddp <- data.frame(
    Size = c(olddata$Size, newdata$Size),
    initial = c(olddata[[init_dens]], newdata[[init_dens]]),
    killed = c(olddata[[response]], rep(NA_integer_, n_new)))

## set up knots
sfun <- function(bs = "cv")  {
    data.frame(Size = smoothCon(s(Size, bs = bs), data = x, absorb.cons = TRUE)[[1]]$knots)
}
k_ar <- sfun()
k_ht <- sfun(bs = "mpd")

ee <- fit$obj$env
pp <- c(split(unname(ee$last.par.best), names(ee$last.par.best)))
ci_level <- 0.95

preds0 <- fit_semimech_fun(data = ddp,
                      response = response,
                      size = ddp[[init_dens]],
                      knots = list(k_ar, k_ht),
                      predict = TRUE,
                      parms = pp)
qq <- qnorm(0.975)
RTMB_sm_pred <- (preds0
    |> filter(nm == "eta")
    |> slice_tail(n = n_new)
    |> transmute(initial = newdata$initial, size = newdata$Size, prop = plogis(value), lwr = plogis(value-qq*sd),
                 upr = plogis(value+qq*sd))
)

odo_sm_aic <- get_info(m_RTMB_sm, newdata = x, init_dens = "initial")
save(m_RTMB_sm, RTMB_sm_pred, odo_sm_aic, file = "odo_semimech.rda")
