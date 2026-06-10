## Load all neccessary libraries
library(gdata)
library(ggplot2)
library(bbmle)
library(R2WinBUGS)
library(plotrix)
library(nlme)
library(MASS)
library(optimx)

do.plots <- FALSE
do.fits <- TRUE
## wd <- "C:/Documents and Settings/anura/Desktop/A callidryas research/Response Surface Experiment/"
wd <- getwd()
datfn <- "McCoy_response_surfaces_Gamboa.csv"

##Import data
## 
x <- read.csv(file.path(wd,datfn))

x$block <- factor(x$block) #Make the block labels into factors

##Visually explore the data
x <- drop.levels(subset(x,cohort=="single"))  # Isolate the data from the single cohort treatments for analysis

bugdat <-  subset(x,predtype=="belo")  #Isolate data for Odonates

bugdat$csize <- bugdat$size - mean(bugdat$size)

## make everything more compact for fitting ...
tmpfun <- function(base,h,size="size",bound=TRUE) {
  basefun <- switch(base,
                    ricker= "1/(c*(size/d)*exp(1-size/d))",
                    powricker="1/(c*(size/d)^gamma*(exp(1-size/d)))",
                    logist= "1/(c/(1+exp((size-d)/b)))",
                    hyper = "1/(c/(1+size/d))",
                    exp = "1/(c*exp(1-size/d))",
                    classic = "1/c")
  hfun <- switch(h,
                 exp= paste("h*exp(hS*",size,")",sep=""),
                 lin= paste("(h+hS*",size,")",sep=""),
                 prop= paste("h*",size,sep=""),
                 ind = "h")
  ffun <- paste("1/(",basefun,"+",hfun,"* initial)",sep="")
  if (bound) {
    ffun <- paste("pmax(eps,pmin(1-eps,",ffun,"))",sep="")
  }
  as.formula(paste("killed~dbinom(prob=",ffun,",size=initial)",sep=""))
}

sv <- list(c=.25,d=3,h=0.05) 
startfun <- function(base,h) {
  c(sv,
    switch(base,
           classic=,ricker=,hyper=,exp=NULL,
           powricker=list(gamma=1),
           logist=list(b=10)),
    switch(h,
           ind=,prop=NULL,
           exp=,lin=list(hS=0)))
}

lfun <- function(base,h,minval=0.002) {
  c(c(c=minval,d=minval,h=minval),
    switch(base,
           classic=,ricker=,hyper=,exp=NULL,
           powricker=c(gamma=-2),
           logist=c(b=0)),
    switch(h,
           ind=,prop=NULL,
           exp=,lin=c(hS=-10)))
}


eps <- 1e-4
## go ahead and fit EVERYTHING
mlist <- c("powricker","exp","ricker","classic","logist","hyper")
## hlist <- c("exp","lin","prop","ind")
hlist <- c("exp","prop","ind") ## only take non-negative models
## olist <- c("BFGS","Nelder-Mead") ## previously: "CG" as option
olist <- c("L-BFGS-B","bobyqa") ## for use with optimizer="optimx"; constraints
slist <- c("size","csize")

k <- 1
res <- list()

## m <- mlist[1]; h <- hlist[1]; o <- olist[2]; s <- slist[1]
## res[[k]] <- mle2(tmpfun(m,h,s),start=startfun(m,h),data=bugdat,
##                     optimizer="optimx",method=o,control=list(maxit=1000),
##                     lower=c(c=0.005,d=0.005,h=0.005,gamma=0.005,hS=-10))

resframe <- expand.grid(sizevar=slist,optmethod=olist,hmodel=hlist,model=mlist)

if (do.fits) {
    for (m in mlist) {
        for (h in hlist) {
            for (o in olist) {
                for (s in slist) {
                    cat(m,h,o,s,"\n")
                    res[[k]] <-  try(
                        mle2(tmpfun(m,h,s),start=startfun(m,h),data=bugdat,
                             optimizer="optimx",method=o,control=list(maxit=1000),
                             lower=lfun(m,h,0))
                    )
                    k <- k+1
                }
            }
        }
    }
    names(res) <- apply(as.matrix(resframe),1,function(x) paste(rev(x),collapse="."))
    save("resframe","res",file="waterbug_fits_1.RData")
    save("resframe","res",file="waterbug_fits_2.RData") ## ?? different answers -- diff optimx version?
    save("resframe","res",file="waterbug_fits_2z.RData") ## zero-boundary case
}
## load("waterbug_fits_1.RData")

## what is this fit from?
load("waterbug_fits_2z.RData")
## some failures
badfit <-  sapply(res,class)=="try-error"
sum(badfit) ## 4 ?  11 ?
resframe[badfit,] ## L-BFGS-B with prop hmodel, csize
res <- res[!badfit] ## drop failures
AICtab(res,sort=TRUE)
rframe <- cbind(resframe[!badfit,],AIC=sapply(res,AIC))

if (do.plots) {
  g1 <- ggplot(rframe,aes(x=hmodel,y=AIC,colour=optmethod,pch=sizevar))+
    geom_point()+facet_grid(.~model)+theme_bw()
  g1

  bestAIC <- min(rframe$AIC)
  g2 <- g1+ylim(bestAIC,bestAIC+20)
  g2

  g2 + geom_hline(yintercept=bestAIC+c(0,2),colour="lightblue",lty=2)
}
rr <- res[rframe$AIC<min(rframe$AIC)+10]
AICtab(rr, sort=TRUE,weights=TRUE)

## source("coefplot_new.R")
## coeftab(res)

## check out best fit
w <- which.min(rframe$AIC)
bestfit <- res[[w]]
rframe[w,]

## everything below changes if we use zero boundary: ignore for now ...
## (i.e., different 'best fit', no longer worrying about ricker/ind ...
##  because probably an artifact of an h=0.002 boundary ...
###
coef(bestfit)  ## ?? minimizes h (on boundary)
## ?? this should be the same as csize/bobyqa/ind/ricker?

fit2 <- res[["ricker.ind.bobyqa.csize"]]
coef(fit2)

## overall inverse-attack rate is (1/(a+hN))
## look at hN and 1/a predictions
with(as.list(coef(bestfit)),curve(1/(c*(x/d)*exp(1-x/d)),from=0,to=50,
                                  ylab="inverse prob [1/a]",ylim=c(0,50),
                                  main="proportional h model"))
## add (hN) prediction for N=1
with(as.list(coef(bestfit)),curve(h*x*1,add=TRUE,col=1,lty=2))
with(as.list(coef(bestfit)),curve(h*x*100,add=TRUE,col=2,lty=2)) ## high-density

x11() ## new window
with(as.list(coef(fit2)),curve(1/(c*(x/d)*exp(1-x/d)),from=0,to=50,
                                  ylab="inverse prob [1/a]",ylim=c(0,50),
                               main = "constant h model"))
## add (hN) prediction for N=1
with(as.list(coef(fit2)),curve(h*1+0*x,add=TRUE,col=1,lty=2))
with(as.list(coef(fit2)),curve(h*100+0*x,add=TRUE,col=2,lty=2)) ## high-density

logLik(fit2)
logLik(bestfit)
## why are these different?

## It should be the case that the h='prop'  and h='ind' models converge
## when h=0 (because one has the term h*initial and the other has the
## term h*csize*initial.  However, I'm getting different answers at
## least in the previous case where the lower parameters were bounded at 0.002
## instead of 0 (a normal precaution in the case where the model is ill-behaved
## outside the boundaries, and where one is using L-BFGS-B, which is particularly
## sensitive to NAs, but maybe not so necessary for this case (first, the model
## isn't necessarily insane for h slightly < 0, secondly, it looks like bobyqa
## might be a bit more robust/respectful of the parameter boundaries than L-BFGS-B
## is -- the latter often violates the boundaries slightly in computing finite-difference
## approximations.

## I'm also getting warnings from optimx that 'parameters or bounds appear to have
## different scalings' -- presumably h is the culprit (small changes away from zero
## might have relatively large effects on the model?), although it would be nice
## if optimx were able to give us some hints about the scaling?  Might be able to
## resolve/fix by adding parscale=c(c=1,d=1,h=0.001) -- but see below

## I'm going to take a stab at getting WinBUGS/JAGS going, instead.  It would
## still be nice to follow this up ...

## now we don't want to bother with this -- we've decided that the
##  ricker fit being best was a bogus consequence of setting the lower
##  bound at h=0.002 rather than h=0 ...
if (FALSE) {
## set things up to fit size-independent model
m <- "ricker"; h <- "ind"; s <- "csize"; o <- "bobyqa"
newfit <- mle2(tmpfun(m,h,s),start=as.list(coef(bestfit)),data=bugdat,
                   optimizer="optimx",method=o,control=list(maxit=1000),
                   lower=lfun(m,h))  ## default lower bound: 0.002

newfit2 <- mle2(tmpfun(m,h,s),start=as.list(coef(bestfit)),data=bugdat,
                   optimizer="optimx",method=o,control=list(maxit=1000),
                   lower=lfun(m,h,0)) ## set the lower bound to 0

## maybe I can get rid of the parameter-scaling warning by explicitly setting parscale
newfit4 <- mle2(tmpfun(m,h,s),start=as.list(coef(bestfit)),data=bugdat,
                   optimizer="optimx",method=o,control=list(maxit=1000,
                                                 parscale=c(c=1,d=1,h=1e-2)),
                   lower=lfun(m,h,0))
## gives error 'pars' not found.  Bug in optimx???

## run proportional model allowing h to go all the way to 0
h <- "prop"
newfit3 <- mle2(tmpfun(m,h,s),start=as.list(coef(bestfit)),data=bugdat,
                   optimizer="optimx",method=o,control=list(maxit=1000),
                   lower=lfun(m,h,0))
logLik(newfit3) ## prop
logLik(newfit2) ## prop, bound at 0
logLik(newfit) ## prop, bound at 0.002
## TO DO: get profile working for 'newfit2' !
##  probably (?) a problem with optimx integration?


## checking models ...
tmpfun("ricker","prop","csize")
tmpfun("ricker","ind","csize")

## urgh.  Work on this later (again!)
## this is a way forward but I haven't finished figuring out what it means
library(emdbook)
z <- calcslice(fit2,bestfit)

plot(z)
abline(h=c(-logLik(fit2),-logLik(bestfit)),col=c(2,4))
z2 <- calcslice(bestfit,fit2)

plot(z2)
abline(h=c(-logLik(fit2),-logLik(bestfit)),col=c(2,4))
}  ## skip evaluation of ricker function with prop/ind etc...
#################
best4 <- rframe$AIC<min(rframe$AIC)+4

rr2 <- res[best4] ## all fits with dAIC<4
## reorder by delta-AIC
rr2 <- rr2[order(sapply(rr2,AIC))]

rframe <- rframe[best4,]
rframe <- rframe[order(rframe$AIC),]

## check sync between rframe and rr2
head(rframe)
head(names(rr2))
## lapply(rr2,coef)
rframe$dAIC <- rframe$AIC-min(rframe$AIC)
rframe$h <- sapply(rr2,function(x)coef(x)["h"])

### compare predictions of ind/powricker with ind/logistic
### how different are they??

tmpfun("logist","ind","size")
tmpfun("powricker","ind","size")
subset(rframe,hmodel=="ind" & model=="powricker") ## check equivalence: L-BFGS-B is actually better

sizecats <- with(bugdat,tapply(size,list(sizeclass),mean))

predframe <- expand.grid(sizeclass=1:5,
                         initial=1:99)
predframe$size <- sizecats[predframe$sizeclass]
predframe$sizeclass <- factor(predframe$sizeclass)
logistpred <- predict(rr2[["logist.ind.L-BFGS-B.size"]],newdata=predframe)
powrickpred <- predict(rr2[["powricker.ind.L-BFGS-B.size"]],newdata=predframe)

pframe2 <- cbind(predframe,logist=logistpred,powricker=powrickpred)
pframe3 <- reshape2::melt(pframe2,id.var=1:3)
pframe3$g <- with(pframe3,interaction(sizeclass,variable))

bugdat$sizeclass <- factor(bugdat$sizeclass)
ggplot(pframe3,aes(x=initial,y=value,colour=sizeclass,group=g,lty=variable))+
  geom_line()+
  geom_point(data=bugdat,aes(y=killed,group=NULL,lty=NULL))

## redo with continuous sizes and discrete initial values

predframe <- expand.grid(size=5:60,
                         initial=c(6,9,18,27,66,81,99))
logistpred <- predict(rr2[["logist.ind.L-BFGS-B.size"]],newdata=predframe)
powrickpred <- predict(rr2[["powricker.ind.L-BFGS-B.size"]],newdata=predframe)

pframe2 <- cbind(predframe,logist=logistpred,powricker=powrickpred)
pframe3 <- reshape2::melt(pframe2,id.var=1:2)
pframe3$g <- with(pframe3,interaction(factor(initial),variable))

bugdat$initial <- factor(bugdat$initial)
if (do.plots) {
ggplot(pframe3,aes(x=size,y=value,colour=factor(initial),group=g,lty=variable))+
  geom_line()+
  geom_point(data=bugdat,aes(y=killed,group=NULL,lty=NULL))+facet_wrap(~initial)
}

  ## new stuff
load("waterbug_fits_2z.RData")  ## load zero-bounded fits
PRzfit <- res[["powricker.ind.L-BFGS-B.size"]]
Lzfit <- res[["logist.ind.L-BFGS-B.size"]]


## do whatever you have to do to load bugdat here ....

m <- "powricker"; h <- "ind"; s <- "size"; o <- "L-BFGS-B"
powrick_ind_new <- mle2(tmpfun(m,h,s),
                        start=as.list(coef(PRzfit)),
                        data=bugdat,
                   optimizer="optimx",method=o,control=list(maxit=1000),
                   lower=lfun(m,h,-0.1))  ## default lower bound: 0.002

m <- "logist"
logist_ind_new <- mle2(tmpfun(m,h,s),start=as.list(coef(Lzfit)),data=bugdat,
                   optimizer="optimx",method=o,control=list(maxit=1000),
                   lower=lfun(m,h,-0.1))  ## default lower bound: 0.002

############
p1 <- profile(powrick_ind_new,which="h")
p2 <- profile(logist_ind_new,which="h") ## uh-oh, found better fit
coef(p2)
logist_ind_new2 <- mle2(tmpfun(m,h,s),start=as.list(coef(p2)),data=bugdat,
                        optimizer="optimx",method=o,
                        control=list(maxit=1000,trace=TRUE),
                        lower=lfun(m,h,-0.1))
## says 'infeasible'? did parameter order get scrambled somehow in bounds?


coef(p2)
p2B <- profile(logist_ind_new2,which="h") ## uh-oh, found better fit

