library(tidyverse); theme_set(theme_bw())
library(plotly)
library(bbmle)
library(emdbook)
library(mgcv)
library(scam)
if (packageVersion("scam") < "1.2.17.9000") stop("need up-to-date/hacked version, see BMB github")
library(RTMB)

source("funs.R")
## rayshader is disappointing/frustrating
## library(rayshader)

## other 3d options? rgl, scatterplot3d, ... ?

datfn <- "McCoy_response_surfaces_Gamboa.csv"

## odonates only
x <- (read.csv(datfn)
    |> transform(block = factor(block))
    |> subset(cohort == "single" & predtype == "odo",
              select = -c(cohort, predtype))
    |> droplevels()
    |> transform(csize = size - mean(size), prop = killed/initial)
)

## ar = f(c(block), size)
## h  = f(h*size)
## prob = f(ar, h)

gg1 <- ggplot(x, aes(x = initial, y = size)) +
    geom_point(aes(colour = killed)) +
    scale_colour_viridis_c()

gg2 <- ggplot(x, aes(x = initial, y = size)) +
    ## stat_sum(aes(colour= killed/initial)) +
    geom_point(aes(colour = killed/initial), size = 4) +
    scale_colour_viridis_c()
## print(gg2)

marker <- list(color = ~prop,
               colorscale = c('#FFE1A1', '#683531'), 
               showscale = TRUE)

seg_data <- function(x, zvar) {
    ## why do we need enquo here??
    zvar <- enquo(zvar)
    xx <- (x
        |> mutate(.id = seq(nrow(x)))
        |> reframe(!!zvar := c(!!zvar, 0), across(-!!zvar), .by = .id)
        |> plotly::group2NA(".id")
    )
    return(xx)
}

## https://stackoverflow.com/questions/72281954/keep-other-columns-when-doing-group-by-summarise-with-dplyr
## https://stackoverflow.com/questions/50012328/r-plotly-showlegend-false-does-not-work
odo_plotly_0 <- (plot_ly(x= ~initial, y = ~size, z = ~prop)
    |> add_markers(data = x, marker = marker, showlegend = FALSE)
    |> add_paths(data = seg_data(x, prop), , showlegend = FALSE)
    |> hide_colorbar()
    |> layout(scene = list(yaxis = list(rangemode = "tozero"),
                           xaxis = list(rangemode = "tozero"),
                           camera = list(eye = list(x = 2.5, y = 2, z = 1)),
                           showlegend=FALSE))
)
if (interactive()) print(odo_plotly_0)

img <- function(obj, width = 1000, height = 1000, trim = TRUE) {
    nm <- paste0(deparse(substitute(obj)), ".png")
    save_image(obj, nm, width = width, height = height)
    if (trim) {
        system(sprintf("convert %s -trim tmp.png; mv tmp.png pix/%s", nm, nm))
    }
}
img(odo_plotly_0)

## https://www.datanovia.com/en/blog/how-to-create-a-ggplot-like-3d-scatter-plot-using-plotly/
## https://community.plotly.com/t/droplines-from-points-in-3d-scatterplot/4113/10

## are these for belo or odo??
L <- load("waterbug_fits_2.RData")
L <- load("waterbug_fits_2z.RData")
aictab <- tibble(resframe, aic = sapply(res, AIC)) |> mutate(across(aic, ~ . - min(.))) |> arrange(aic)
    

## power-ricker or  ricker attack rate, proportional handling time
m_mle2_rickerprop <- mle2(killed ~ dbinom(prob = 1/(1/(c*size/d*exp(1-size/d)) + h*size*initial),
                             size = initial),
             start = list(c=1,d=20,h=20),
             control = list(parscale= c(c=1,d=20,h=20)),
             data = x,
             method = "Nelder-Mead")

long_fmt <- function(cc, nm = c("size", "initial", "prop")) {
    tibble(x = rep(cc$x, ncol(cc$z)),
           y = rep(cc$y, each = nrow(cc$z)),
           z = c(cc$z)) |>
        setNames(nm)
}
    
cc <- curve3d(1/(1/(c*size/d*exp(1-size/d)) + h*size*initial),
              xlim = c(0, 60),
              ylim = c(0, 100),
##            1/(1/(c*size/d*exp(1-size/d)) + h*size*initial)
        data = as.list(coef(m_mle2_rickerprop)),
        varnames = c("size", "initial"),
        sys3d = "image")


## https://stackoverflow.com/questions/34178381/how-to-specify-camera-perspective-of-3d-plotly-chart-in-r
odo_plotly_param <- (odo_plotly_0
    |> add_trace(type =  "mesh3d", data = long_fmt(cc), opacity = 0.4)
)
img(odo_plotly_param)

m_gam_te <- gam(cbind(killed, initial-killed) ~ te(size, initial),
            data = x, family = binomial)

odo_gam_pred <- expand.grid(size = 0:60, initial = 0:100)
odo_gam_pred$prop <- predict(m_gam_te, newdata = odo_gam_pred, type = "response")

odo_plotly_gam <- (odo_plotly_0
    |> add_trace(type =  "mesh3d", data = odo_gam_pred, opacity = 0.4)
)
img(odo_plotly_gam)

xx <- x |> select(killed, size, initial) |> expand_bern(response = "killed", size = x$initial)
## tesmd2, tesmd1 = smooth monotone decreasing in var1/2, no constraint otherwise
## fit3 <- scam(killed ~ s(initial, size, bs = "tesmd2"), data = xx, family = binomial)
## fit3 <- scam(killed ~ s(initial, size, bs = "tesmd1"), data = xx, family = binomial)
## fails for belo .... (inner loop 3, can't correct step size)

## https://plotly.com/r/static-image-export/

## decreasing wrt var 1, convex wrt var 2
m_scam_tedecv <- scam(killed ~ s(initial, size, bs = "tedecv"), data = xx, family = binomial)

scam_pred <- odo_gam_pred
scam_pred$prop <- predict(m_scam_tedecv, newdata = scam_pred, type = "response")
odo_plotly_scam <- odo_plotly_0 |> add_trace(type =  "mesh3d", data = scam_pred, opacity = 0.4)
img(odo_plotly_scam)

ss <- s(initial, size, bs = "tedecv")
x2 <- x |> rename(Size = "size") ## hack, 'size' is confounded
m_RTMB_tedecv <- fit_mpd_fun(data = x2[c("killed", "Size", "initial")],
            response = "killed",
            xvar = c("Size", "initial"),
            form = s(size, initial, bs = "tesmd2"),
            size = x2$initial,
            family = "binomial",
            inner.control = list(smartsearch=FALSE, maxit =1),
            opt = "BFGS",
            start = list(b0 = -2, log_smD = 2),
            ## works with random = NULL; start from better values?
            ## 'initial value in vmmin is not finite' ...
            random = NULL)
m_RTMB_tedecv$fit

c(nll_gam = -sum(dbinom(x2$killed, size = x2$initial, prob = predict(m_gam_te, type = "response"), log = TRUE)),
  nll_RTMB = -sum(dbinom(x2$killed, size = x2$initial, prob = m_RTMB_tedecv$mu, log = TRUE)))

get_info(m_scam_tedecv, newdata = x, init_dens = "initial")
get_info(m_gam_te, newdata = x, init_dens = "initial")
get_info(m_RTMB_tedecv, newdata = x, init_dens = "initial")

## adapt rf_predfun guts rather than trying to make everything universal/back-compatible
newdata <- odo_gam_pred
olddata <- x2
init_dens <- "initial"
response <- "killed"
fit <- m_RTMB_tedecv
n_new <- nrow(newdata)

ddp <- data.frame(
    Size = c(olddata$Size, newdata$size),  ## change cap again
    initial = c(olddata[[init_dens]], newdata[[init_dens]]),
    killed = c(olddata[[response]], rep(NA_integer_, n_new)))
k <- data.frame(initial = smoothCon(s(Size, initial, bs="tesmd2"), data = olddata, absorb.cons = TRUE)[[1]]$knots)
ee <- fit$obj$env
pp <- c(split(unname(ee$last.par.best), names(ee$last.par.best)))
ci_level <- 0.95

preds0 <- fit_mpd_fun(data = ddp, response = response,
                      size = ddp[[init_dens]],
                      xvar = c("Size", "initial"),
                      form = s(x, y, bs = "tesmd2"), 
                      family = "binomial",
                      random = NULL,
                      knots = k,
                      predict = TRUE,
                      parms = pp)
qq <- qnorm((1+ci_level)/2)
preds_RTMB <- (preds0
    |> filter(nm == "eta")
    |> slice_tail(n = n_new)
    |> transmute(initial = newdata$initial, size = newdata$size, prop = plogis(value), lwr = plogis(value-qq*sd),
                 upr = plogis(value+qq*sd))
)
odo_plotly_RTMB_tecdv <- odo_plotly_0 |> add_trace(type =  "mesh3d", data = preds_RTMB, opacity = 0.4)
print(odo_plotly_RTMB_tecdv)
img(odo_plotly_RTMB_tecdv)

## include m_RTMB_tedecv ???
odo_aictab <- (map_dfr(tibble::lst(m_scam_tedecv, m_gam_te, m_mle2_rickerprop),
                      \(m) get_info(m, newdata = x, init_dens = "initial"),
                      .id = "model")
    |> arrange(AIC)
    |> mutate(across(c(AIC, nll), ~ . - min(., na.rm = TRUE)))
    |> rename_with(\(x) paste0("Δ", x), c(AIC, nll))
    |> mutate(across(model, \(x) gsub("_", "/", gsub("^m_", "", x))))
)

save("odo_aictab", file = "odo_stuff.rda")

load("odo_semimech.rda")
odo_plotly_RTMB_sm <- odo_plotly_0 |> add_trace(type =  "mesh3d", data = RTMB_sm_pred, opacity = 0.4)
print(odo_plotly_RTMB_sm)
img(odo_plotly_RTMB_sm)

## original waterbug JAGS model for odonates (power-Ricker + prop, random effects of block)
## 
## for (i in 1:N) {
##     ## ar[i] <- cvec[block[i]]*pow(size[i]/d,gamma)*exp(1-size[i]/d)
##     ar[i] <- cvec[block[i]]*size[i]/d*exp(1-size[i]/d)
##     ## hfun[i] <- h*exp(hS*size[i])
##     hfun[i] <- h*csize[i]  ## prop.
##     prob[i] <- 1/(1/ar[i]+hfun[i]*initial[i])
##     killed[i] ~ dbin(prob[i],initial[i])
##   }
##   for (i in 1:nblock) {
##     cvec0[i] ~ dnorm(0,tau.c)
##     cvec[i] <- c*exp(cvec0[i])
##     killed.rep[i]~ dbin(prob[i],initial[i])
##   }
## ## priors
##   tau.c <- pow(sd.c,-2)
##   sd.c ~ dunif(0,1)
##   d ~ dlnorm(0,0.01)
##   gamma ~ dlnorm(0,0.01)
##   h ~ dunif(0,10) ## changed to positive uniform
##   ## hS ~ dunif(-1,1) 
##   c ~ dlnorm(0,0.01)

