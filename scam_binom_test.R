## make sure we have the latest/bug-fixed versions
remotes::install_github("bbolker/reformulas")
remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

library(RTMB)
library(glmmTMB)
library(scam)
source("funs.R")

if (!interactive()) pdf("scam_binom_test.pdf")

set.seed(101)
dd <- data.frame(x=seq(-5, 5, length = 101))
dd <- within(dd, {
             p <- plogis(1 + x - x^3/4)
             y0 <- rbinom(nrow(dd), size = 1, prob = p)
             y1 <- rbinom(nrow(dd), size = 20, prob = p)
             })

par(las = 1, bty = "l") ## cosmetic
colvec <- c(2,4:6, 8:9)

## ONLY works without Laplace approx ... maybe OK to proceed for now?
m_bern_RTMB <- fit_mpd_fun(dd, size = rep(1, nrow(dd)), family = "binomial", response = "y0",
                           random = NULL,
                           ## random = "b1",
                           parms = list(b0 = 0, log_smSD = 3, b1 = rep(0, 9)))

## all approaches behave similarly for Bernoulli data
## scam only does GCV, glmmTMB only does ML/REML
m_bern_gam_gcv <- gam(y0 ~ s(x, bs = "tp"), family = binomial, data = dd,
                      method = "GCV.Cp")
m_bern_gam_reml <- gam(y0 ~ s(x, bs = "tp"), family = binomial, data = dd,
                       method = "REML")
m_bern_scam <- scam(y0 ~ s(x, bs = "tp"), family = binomial, data = dd)
m_bern_glmmTMB <- glmmTMB(y0 ~ s(x, bs = "tp"), family = binomial, data = dd,
                          REML = TRUE)

nm0 <- ls(pattern = "m_bern_.*")
nm <- nm0 |> gsub(pattern = "^m_", replacement = "")
predmat_bern <- mget(nm0) |> lapply(predict) |> do.call(what=cbind)
colnames(predmat_bern) <- nm

matplot(dd$x, predmat_bern, lty = 1:4, col = colvec)
legend(lty = 1, lwd = 2,
       col = colvec[1:4],
       x = 1.5, y = 8,
       legend = c("gam/GCV", "gam/REML", "scam/GCV", "glmmTMB/REML"))

lines(dd$x, qlogis(dd$p), lwd = 2)
rug(dd$x[dd$y0==0], ticksize = 0.1, side = 1)
rug(dd$x[dd$y0==1], ticksize = 0.1, side = 3)
## results are heavily smoothed relative to true prob, but not surprising given
## v. low-resolution data

m_binom_gam_gcv <- gam(cbind(y1, 20-y1) ~ s(x, bs = "tp"), family = binomial, data = dd,
                      method = "GCV.Cp")
m_binom_gam_reml <- gam(cbind(y1, 20-y1) ~ s(x, bs = "tp"), family = binomial, data = dd,
                       method = "REML")
m_binom_scam_2col <- scam(cbind(y1, 20-y1) ~ s(x, bs = "tp"), family = binomial, data = dd)
## NOTE warnings() about non-integer # successes (due to weights issues)
m_binom_scam_wts <- scam(y1/20 ~ s(x, bs = "tp"), weights = rep(20, nrow(dd)),
                                                             family = binomial, data = dd)
m_binom_glmmTMB <- glmmTMB(cbind(y1, 20-y1) ~ s(x, bs = "tp"), family = binomial, data = dd,
                           REML = TRUE)

dd_expand <- expand_bern(dd)
m_binom_scam_expand <- scam(y1 ~ s(x, bs = "tp"),
                            family = binomial, data = dd_expand)

nm0 <- ls(pattern = "m_binom_.*")
nm <- nm0 |> gsub(pattern = "^m_", replacement = "")

## RTMB silently ignores newdata. no-op for the rest, only relevant for expanded data
predmat_binom <- (mget(nm0)
    |> lapply(predict, newdata = dd)
    |> as.data.frame()
    |> setNames(nm)
)

matplot(dd$x, predmat_binom, lty = 1:4, col = colvec, ylim = c(-25,25))
lines(dd$x, qlogis(dd$p), lwd = 2)
## squeeze obs probs in slightly so we can look at logit scale
points(dd$x, qlogis((dd$y1+0.25)/20.5))
legend("topright", lty = 1, lwd = 2,
       col = colvec,
       legend = gsub("_", "/", gsub("binom_", "", nm))
       )

## note, RTMB is mpd, not tp; FIXME ...

m_binom_RTMB_mpd <- fit_mpd_fun(dd, size = rep(20, nrow(dd)), family = "binomial", response = "y1",
                            random = "b1",
                            parms = list(b0 = 0, log_smSD = 3, b1 = rep(0, 9)))
m_binom_scam_expand_mpd <- scam(y1 ~ s(x, bs = "mpd"),
                            family = binomial, data = dd_expand)

nm0 <- ls(pattern = "m_.*_mpd")
nm <- nm0 |> gsub(pattern = "^m_", replacement = "")

## RTMB silently ignores newdata. no-op for the rest, only relevant for expanded data
predmat_mpd <- (mget(nm0)
    |> lapply(predict, newdata = dd)
    |> as.data.frame()
    |> setNames(nm)
)

matplot(dd$x, predmat_mpd, lty = 1:2, type = "l", col = colvec, ylim = c(-25,25))
lines(dd$x, qlogis(dd$p), lwd = 2)
## squeeze obs probs in slightly so we can look at logit scale
points(dd$x, qlogis((dd$y1+0.25)/20.5))
legend("topright", lty = 1, lwd = 2,
       col = colvec,
       legend = gsub("_", "/", gsub("binom_", "", nm))
       )

## close, not identical. RTMB actually looks slightly more sensible (ML vs GCV difference?)

if (!interactive()) dev.off()
