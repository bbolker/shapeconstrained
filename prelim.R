library(scam)
library(glmmTMB)
library(marginaleffects)
library(ggplot2); theme_set(theme_bw())
library(purrr)
library(dplyr)
## what for tidy predictions? ggpredict, marginaleffects, emmeans?
source("funs.R")

## 1. simulate data from Holling type-2, type-3:
set.seed(101)

data("ReedfrogFuncresp", package = "emdbook")
dd <- ReedfrogFuncresp
ddx <- expand_bern(dd, response = "Killed", size = "Initial")
## apropos("smooth.construct")
## ?smooth.construct.miso.smooth.spec
## scam smooth codes:
##  m = monotonic
##  p = p-spline
##  i/d = increasing/decreasing
##  cv = concavity
## te = tensor
## d = double
## de = 'decreasing' (why not md?)

## linear fits
scam_rf_mpd_gauss <- scam(Killed/Initial ~ s(Initial, bs = "mpd"), data = dd)
gam_rf_tp_gauss <- gam(Killed/Initial ~ s(Initial, bs = "tp", k = 8), data = dd)
gam_rf_tp_binom <- gam(cbind(Killed, Initial -Killed) ~ s(Initial, bs = "tp",
                                                          k = 8), data = dd,
                       family = binomial)
glm_rf_holling_binom <- glm(cbind(Killed, Initial -Killed) ~ I(1/Initial),
                            data = dd, family = binomial)


## can't do REML selection ...
scam_rf_mpd_binom <- scam(Killed ~ s(Initial, bs = "mpd"), data = ddx, family = binomial)

predfun <- function(m) {
    nd <- data.frame(Initial = 1:100)
    if (!inherits(m, "glmmTMB")) {
        c(marginaleffects::predictions(m, newdata = nd) |> as_tibble())
    } else stop("can't do glmmTMB preds yet")
}

models <- ls(pattern="^(scam|gam|glm)_rf")
mod_list <- mget(models) |> setNames(models)
preds <- (mod_list
    |> map_dfr(predfun, .id = "model")
    |> select(model, Initial, prob = estimate, lwr = conf.low, upr = conf.high)
)
    

ggplot(preds, aes(Initial, prob)) +
    geom_line(aes(colour = model)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), colour = NA, alpha = 0.5) +
    expand_limits(y=0) +
    geom_point(data=dd, aes(y = Killed/Initial, size = Killed), alpha = 0.5) +
    facet_wrap(~model)
## CIs are not monotonic??


scam_pos <- match("package:scam", search())
aa <- apropos("smooth.construct", where = TRUE)
scam_smooths <- unname(aa[names(aa) == scam_pos]) |>
    gsub(pattern = "[.]?smooth\\.(construct|spec)[.]?", replacement = "")
