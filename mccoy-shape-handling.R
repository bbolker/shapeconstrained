## get data
x <- read.csv("Response surfaces data_Gamboa _data.csv")
library(gdata)
x$block <- factor(x$block)
x <- drop.levels(subset(x,cohort=="single"))

sv <- list(c=0.04,d=40,h=0.025)
## crank down c a little bit to get all initial mort < 1
ododat <-  subset(x,predtype=="odo")

library(bbmle)
## power-Ricker
g3_powricker = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d)^gamma*(exp(1-size/d)))+
                       h*initial),
  size=initial),
  start=c(sv,list(gamma=1)),
  data=ododat)
eps <- 1e-4
g3_powricker_pen <- 
  mle2(killed~dbinom(prob=pmax(eps,
                       pmin(1-eps,
                            1/(1/(c*(size/d)^gamma*(exp(1-size/d)))+
                            h*initial))),
                     size=initial),
       start=c(sv,list(gamma=1)),
       data=ododat)
all.equal(coef(g3_powricker),coef(g3_powricker_pen)) ## ended up in the same place

ododat$csize <- ododat$size - mean(ododat$size)

## Hessian is ending up exactly singular; wants to make gamma arbitrarily
##  large, h arbitrarily small -- i.e. linear fits
sv2 <- as.list(coef(g3_powricker))
sv2$h <- NULL
sv2$logh <- -2

g3_powricker_h1 = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d)^gamma*(exp(1-size/d)))+
                       exp(logh)*size*initial),
                     size=initial),
       method="Nelder-Mead",
       start=sv2,
       trace=TRUE,
       skip.hessian=TRUE,
       control=list(parscale=abs(unlist(sv2))),
       data=ododat)

library(emdbook)
G <- as.list(coef(g3_powricker_h1)) ##kluge

curve3d(initial*(1/(1/(G$c*(size/G$d)^G$gamma*(exp(1-size/G$d)))+
                    exp(G$logh)*size*initial)),
        varnames=c("initial","size"),
        sys3d="rgl",col="gray")
grid3d(c("x+","y+","z"))

## what is the shape of this surface, and how does it depend on h?
## (this is where a 4D 
