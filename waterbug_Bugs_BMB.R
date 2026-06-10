##Load all neccessary libraries
library(gdata)
library(ggplot2)
library(bbmle)
library(R2jags)
library(plotrix)
library(nlme)
library(MASS)
library(coda)

##Import data
x <- read.csv("McCoy_response_surfaces_Gamboa.csv")
x$block <- factor(x$block) #Make the block labels into factors
bugdat <-  subset(x,predtype=="belo")  #Isolate data for Odonates
## str(bugdat)

L <- load("waterbug_fits_1.RData")  ## gets res, resframe
badfit <-  sapply(res,class)=="try-error"
sum(badfit) ## 4
resframe[badfit,] ## L-BFGS-B with prop hmodel, csize
res <- res[!badfit] ## drop failures
rframe <- cbind(resframe[!badfit,],AIC=sapply(res,AIC))


bestfit <- res[[which.min(rframe$AIC)]]
rframe[which.min(rframe$AIC),]
## Prop/Ricker model, using centered size (fitted with bobyqa)
## I am just going to go ahead and ASSUME that this is the best
##   model, ignoring any warning bells that are going off because
##   of the issues commented on in the companion script (waterbug_BMB.R) ...

xstart <- c(as.list(coef(bestfit)),list(sd.c=0.1))
#xstart$h=.0001
set.seed(1001)
##Create a list of starting based on xstart above
xstart2 = list(
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)),
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)),
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)))
## BMB: might be able to use the (poorly documented) perturb.params()
##   function from emdbook here, although this is just fine
## Now create a bugs data file for Winbugs program
xx <- bugdat
xx$csize <- xx$size-mean(xx$size)
xx$block <- as.numeric(xx$block)
xx$killed.p <- xx$killed/xx$P

dat <- c(as.list(subset(xx,select=c(killed,size,csize,initial,block))),
         N=nrow(bugdat), nblock=length(levels(bugdat$block)))

## I would love to be able to automatically insert the appropriate
##   formula into the bugsmodel ...
bugsmodel <- function() {
  for (i in 1:N) {
    ## ar[i] <- cvec[block[i]]*pow(size[i]/d,gamma)*exp(1-size[i]/d)
    ar[i] <- cvec[block[i]]*size[i]/d*exp(1-size[i]/d)
    ## hfun[i] <- h*exp(hS*size[i])
    hfun[i] <- h*csize[i]  ## prop.
    prob[i] <- 1/(1/ar[i]+hfun[i]*initial[i])
    killed[i] ~ dbin(prob[i],initial[i])
  }
  for (i in 1:nblock) {
    cvec0[i] ~ dnorm(0,tau.c)
    cvec[i] <- c*exp(cvec0[i])
    killed.rep[i]~ dbin(prob[i],initial[i])
  }
## priors
  tau.c <- pow(sd.c,-2)
  sd.c ~ dunif(0,1)
  d ~ dlnorm(0,0.01)
  gamma ~ dlnorm(0,0.01)
  h ~ dunif(0,10) ## changed to positive uniform
  ## hS ~ dunif(-1,1) 
  c ~ dlnorm(0,0.01)
}
body(bugsmodel)[[2]][[4]][[2]]
base <- "ricker"
bugsfun <- paste("ar <-",switch(base,
                                ricker=     "cvec[block[i]]*size[i]/d*exp(1-size[i]/d)",
                                powricker=  "cvec[block[i]]*pow(size[i]/d,gamma)*exp(1-size/d)",
                                logist=     "cvec[block[i]]/(1+exp((size[i]-d)/b))",
                                hyper =     "cvec[block[i]]/(1+size[i]/d)",
                                exp =       "cvec[block[i]]*exp(1-size[i]/d)",
                                classic =   "cvec[block[i]]"))
## don't know how to do this
tmpf <- bugsmodel
## body(tmpf)[[2]][[4]][[2]] <- parse(text=bugsfun)
##parse(text=bugsfun)

fn="ricker_prop_bugs.RData"

bfun <- jags
if (!file.exists(fn)) {
  t.wbugs <- system.time(b <- bfun(data=dat,
    inits=c(list(xstart),xstart2),
    n.chains=4,
    n.iter=80000,
    parameters=c(names(xstart),"cvec","killed.rep"),
    model.file=bugsmodel,
    working.directory=getwd()))
  ## dropped debug=TRUE, unknown to jags
                                  
  save("b","t.wbugs",file=fn)
  } else load(fn)

fn <- "" ## prevent overwriting later on


plot(b$BUGSoutput)
b2 <- as.mcmc.list(b$BUGSoutput)
devpos <- which(colnames(b2[[1]])=="deviance")
krpos <- grep("killed",colnames(b2[[1]]))
b3 <- b2[,-c(devpos,krpos)]
xyplot(b3,layout=c(3,4),asp="fill")
densityplot(b3,layout=c(3,4),asp="fill")
gelman.diag(b3)

source("coefplot_new.R")
print(coeftab(b3),digits=3)
