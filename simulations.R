library(glmmTMB)
library(RTMB)
library(ggplot2); theme_set(theme_bw())
library(purrr)
library(dplyr)
source("funs.R")

## binom mpd example works, reedfrog fails.
## investigate.

## when does Laplace approx work?
## when does RTMB give similar answers to scam?

hpars <- list(a=0.5, h = 1/100)
##' @param n total number of observations
##' @param xmin minimum starting density
##' @param xmax maximum starting density
##' @param reps samples per initial density (n/reps must be integer)
##' @param pars parameters (list: a, h for Holling, b0, b1 for spline)
##' @param random-number seed
simfun <- function(n = 100, xmin=5, xmax=100, reps = 1,
                   pars = hpars, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (n %% reps != 0) warning("non-integer n/reps")
    n_dens <- n/reps
    dd <- expand.grid(Initial = round(seq(xmin, xmax, length.out = n_dens)),
                      rep = seq(reps))
    if ("a" %in% names(pars)) {
        prob <- with(pars, a/(1+a*h*dd$Initial))
    } else {
        sm1 <- smoothCon(s(Initial, bs = "mpd"), data = dd, absorb.cons = TRUE)[[1]]
        if (is.null(pars$b1)) b1 <- rnorm(ncol(sm1$X))
        prob <- with(pars,
                     plogis(b0 + sm1$X %*% exp(b1)))
    }
    dd$Killed <- rbinom(nrow(dd), prob = prob, size = dd$Initial)
    return(dd)
}

set.seed(101)
sims1 <- purrr::map_dfr(1:16,
               \(x) simfun(n = 50),
               .id = "sim")

ggplot(sims1, aes(Initial, Killed/Initial)) +
    geom_point(aes(size = Killed)) + facet_wrap(~sim) + geom_smooth() +
    geom_function(fun = function(x) 0.5/(1+0.5*(1/100)*x), col = "red")

simfun(n=50, pars = list(b0 = 0))
sims2 <- purrr::map_dfr(1:16,
               \(x) simfun(n = 50, pars = list(b0 = 0)),
               .id = "sim")
ggplot(sims2, aes(Initial, Killed/Initial)) +
    geom_point(aes(size = Killed)) + facet_wrap(~sim) + geom_smooth()

## now try to fit Holling to simulated data ...

combfun <- function(..., pars = hpars, random = "b1") {
    dd <- simfun(..., pars = pars)
    m_RTMB_mpd <- fit_mpd_fun(data = dd, response = "Killed",
                              size = dd$Initial, xvar = "Initial",
                              family = "binomial", random = random)
    ddx <- expand_bern(dd, response = "Killed", size = "Initial")
    m_scam_mpd <- m_scam_mpd <- scam(Killed ~ s(Initial, bs = "mpd"), data = ddx,
                                     family = binomial)
    ddp0 <- data.frame(Initial = dd$Initial)
    pfun <- function(m) (c(marginaleffects::predictions(m, newdata = ddp0))
        |> as.data.frame()
        |> as_tibble()
        |> dplyr::select(prob = estimate, lwr = conf.low, upr = conf.high)
        |> mutate(Initial = ddp0$Initial, .before = 1)
    )

    p2 <- pfun(m_scam_mpd)
    true <- if (!is.null(pars$a)) {
                with(pars, a/(1+a*h*dd$Initial))
            } else NULL
            
    tibble(dd,
           scam_prob = p2$prob, scam_lwr = p2$lwr, scam_upr = p2$upr,
           true = true,
           RTMB_prob = drop(m_RTMB_mpd$mu))
}

set.seed(101)
cc <- combfun(random = NULL)
gg0 <- ggplot(cc, aes(x=Initial)) +
    geom_point(aes(y= Killed/Initial)) +
    geom_line(aes(y=scam_prob), colour = "red") +
    geom_line(aes(y=true), colour = "blue") +
    geom_line(aes(y=RTMB_prob), colour = "purple")
print(gg0)

set.seed(101)
sim2 <- purrr::map_dfr(1:16, \(x) combfun(random = NULL),
               .id = "sim")

gg0 %+% sim2 + facet_wrap(~sim)


sim3 <- purrr::map_dfr(1:16, \(x) combfun(random = NULL, pars = list(b0 = 0)),
               .id = "sim")

ggplot(sim3, aes(x=Initial)) +
    geom_point(aes(y= Killed/Initial)) +
    geom_line(aes(y=scam_prob), colour = "red") +
    geom_line(aes(y=RTMB_prob), colour = "purple") +
    facet_wrap(~sim)

