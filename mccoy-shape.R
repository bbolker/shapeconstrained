###################################################
### chunk number 1: double eval=FALSE
###################################################
## #line 19 "mccoy-shape.Rnw"
## options(SweaveHooks=list( fig=function() par(mfrow=c(1,2)) ))


###################################################
### chunk number 2: single
###################################################
#line 23 "mccoy-shape.Rnw"
options(SweaveHooks=list( fig=function() par(mar=c(5,13,4,11)+0.1)) )


###################################################
### chunk number 3: getdata
###################################################
#line 31 "mccoy-shape.Rnw"
x <- read.csv("Response surfaces data_Gamboa _data.csv")
library(gdata)
x$block <- factor(x$block)
x <- drop.levels(subset(x,cohort=="single"))


###################################################
### chunk number 4: design
###################################################
#line 38 "mccoy-shape.Rnw"
with(x,table(sizeclass,initial,predtype))


###################################################
### chunk number 5: loadgg
###################################################
#line 42 "mccoy-shape.Rnw"
library(ggplot2)


###################################################
### chunk number 6: plot1_dens
###################################################
#line 52 "mccoy-shape.Rnw"
#line 19 "mccoy-shape.Rnw#from line#52#"
options(SweaveHooks=list( fig=function() par(mfrow=c(1,2)) ))
#line 53 "mccoy-shape.Rnw"
op=par(xpd=NA)
ss = scale_shape_manual(name="block",values=1:7)
lfun = function(variable,value) {
  ifelse(value=="belo","Water bug","Dragonfly")
}
## testing
p1 = ggplot(x,aes(x=initial,y=killed,group=factor(sizeclass),
             fill=factor(sizeclass),
             colour=factor(sizeclass),
             pch=block))+ss+
  labs(colour="size class",fill="size class")
print(p1 + facet_grid(~predtype,labeller=lfun)+
  stat_sum(aes(size=..n..))+geom_smooth())
par(op)


###################################################
### chunk number 7: plot1_size
###################################################
#line 79 "mccoy-shape.Rnw"
p2 = ggplot(x,aes(x=size,y=killed,group=factor(initial),
  fill=factor(initial),colour=factor(initial),
  pch=block))+labs(colour="initial density",fill="initial density")
insubdat <- subset(x,!factor(initial) %in% c("9","27","81"))
print(p2+geom_point()+facet_grid(~predtype,labeller=lfun)+
      geom_smooth(data=insubdat)+
      coord_cartesian(ylim=c(-5,52)))


###################################################
### chunk number 8: 
###################################################
#line 122 "mccoy-shape.Rnw"
library(bbmle)


###################################################
### chunk number 9: modelfits
###################################################
#line 126 "mccoy-shape.Rnw"
sv <- list(c=0.04,d=40,h=0.025)
## crank down c a little bit to get all initial mort < 1
ododat =  subset(x,predtype=="odo")
g3_ricker = mle2(killed~dbinom(prob=1/(1/(c*(size/d)*exp(-size/d))+h*initial),
  size=initial),
  start=sv,
  data=ododat)
g3_powricker = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d)^gamma*(exp(1-size/d)))+
                       h*initial),
  size=initial),
  start=c(sv,list(gamma=1)),
  data=ododat)
eps = 1e-4
g3_powricker_pen <- 
  mle2(killed~dbinom(prob=pmax(eps,
                       pmin(1-eps,
                            1/(1/(c*(size/d)^gamma*(exp(1-size/d)))+
                            h*initial))),
                     size=initial),
       start=c(sv,list(gamma=1)),
       data=ododat)
g3_genricker <- 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d*exp(1-size/d))^gamma)+
                       h*initial),
  size=initial),
  start=c(sv,list(gamma=1)),
  data=ododat,method="Nelder-Mead",control=list(maxit=1000))

g3_genricker_hS1 = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d*exp(1-size/d))^gamma)+
                       h*size*initial),
  size=initial),
  start=c(sv,list(gamma=1)),
  data=ododat,method="Nelder-Mead",control=list(maxit=1000))
ododat$csize <- ododat$size-mean(ododat$size)
g3_genricker_hS2 = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d*exp(1-size/d))^gamma)+
                       (h+hS*csize)*initial),
  size=initial),
  start=c(sv,list(gamma=1,hS=0)),
  data=ododat,method="Nelder-Mead",control=list(maxit=1000))
g3_logist = mle2(killed~dbinom(prob=1/(1/(c/(1+exp((size-d)/b)))+
                                 h*initial),
  size=initial),
  start=c(sv,list(b=10)),
  data=ododat)
g3_hyper = mle2(killed~dbinom(prob=1/(1/(c/(1+size/d))+
                                 h*initial),
  size=initial),
  start=sv,
  data=ododat)
g3_hyperq = mle2(killed~dbinom(prob=1/(1/(c/(1+(size/d)+phi*(size/d)^2))+
                                 h*initial),
  size=initial),
  ## start=c(as.list(coef(g3_hyper)),list(phi=0)),
  start=list(c=1/5.65,d=(5.65/-0.656),h=0.00698,phi=0.02536*(5.65/(0.656^2))),
  ## start = list(c=4,d=-0.03,phi=2,h=0.013),
  ## method="SANN",control=list(maxit=10000),
  data=ododat)
g3_hyperq_bad = mle2(killed~dbinom(prob=1/(1/(c/(1+(size/d)+phi*(size/d)^2))+
                                 h*initial),
  size=initial),
  start=c(as.list(coef(g3_hyper)),list(phi=0)),
  data=ododat)
g3_exp = mle2(killed~dbinom(prob=1/(1/(c*exp(-size/d))+
                                 h*initial),
  size=initial),
  start=sv,
  data=ododat)
pred.p = predict(g3_powricker)/ododat$initial
resids = with(ododat,
  (predict(g3_powricker)-killed)/sqrt(initial*pred.p*(1-pred.p)))
xr = data.frame(ododat,resid=resids)
disp = sum(resids^2)/(nrow(ododat)-length(coef(g3_powricker)))
(q0=ICtab(g3_ricker,g3_powricker,g3_genricker,
   g3_logist,g3_hyper,g3_hyperq,g3_exp,
   g3_hyperq_bad,
   type="qAIC",sort=TRUE,delta=TRUE,dispersion=disp,weights=TRUE))
modlist = list(powricker=g3_powricker,
  genricker=g3_genricker,
  logist=g3_logist,
  hyperq=g3_hyperq,
  ricker=g3_ricker,
  exp=g3_exp,
  hyper=g3_hyper)
## reorder list to match ICtab
modlist <- modlist[gsub("g3_","",attr(q0,"row.names"))]
## get rid of anything that doesn't match (??)
modlist <- modlist[!is.na(names(modlist))]


###################################################
### chunk number 10: 
###################################################
#line 218 "mccoy-shape.Rnw"
predframe = expand.grid(size=seq(5,60,by=0.5),
  initial=sort(unique(x$initial)))
pfun = function(fit,name) {
  data.frame(predframe,
             killed=predict(fit,newdata=predframe),
             model=factor(name))
}
mframe= do.call(rbind,mapply(pfun,modlist,names(modlist),SIMPLIFY=FALSE))
odo_mframe = data.frame(mframe,predtype="odo")


###################################################
### chunk number 11: plot_models1
###################################################
#line 230 "mccoy-shape.Rnw"
print(ggplot(subset(odo_mframe,model!="hyper"),
             aes(x=size,y=killed,colour=factor(initial),
                 group=factor(initial)))+
      geom_line()+facet_wrap(~model)+geom_point(data=x)+
      labs(colour="initial density"))


###################################################
### chunk number 12: plot_resids2
###################################################
#line 244 "mccoy-shape.Rnw"
print(ggplot(xr,
       aes(x=block,y=resid))+geom_boxplot())
anova(lm(resid~block,data=xr))


###################################################
### chunk number 13: bugmodelfits
###################################################
#line 252 "mccoy-shape.Rnw"
sv = list(c=0.04,d=40,h=0.025)
## crank down c a little bit to get all initial mort < 1
bugdat =  subset(x,predtype=="belo")
g4_ricker = mle2(killed~dbinom(prob=1/(1/(c*(size/d)*exp(-size/d))+h*initial),
  size=initial),
  start=sv,
  data=bugdat)
g4_powricker = mle2(killed~dbinom(prob=1/(1/(c*(size/d)^gamma*exp(1-size/d))+
                                    h*initial),
  size=initial),
  start=c(sv,list(gamma=1)),
  data=bugdat)
g4_genricker = mle2(killed~dbinom(prob=1/(1/(c*(size/d*exp(1-size/d))^gamma)+
                                    h*initial),
  size=initial),
  start=c(sv,list(gamma=1)),
  data=bugdat,method="Nelder-Mead",control=list(maxit=1000))
g4_logist = mle2(killed~dbinom(prob=1/(1/(c/(1+exp((size-d)/b)))+
                                 h*initial),
  size=initial),
  start=c(sv,list(b=10)),
  data=bugdat)
g4_hyper = mle2(killed~dbinom(prob=1/(1/(c/(1+size/d))+
                                 h*initial),
  size=initial),
  start=sv,
  data=bugdat)
g4_hyperq = mle2(killed~dbinom(prob=1/(1/(c/(1+(size/d)+phi*(size/d)^2))+
                                 h*initial),
  size=initial),
  ## start=c(as.list(coef(g4_hyper)),list(phi=0)),
  start=list(c=1/5.65,d=(5.65/-0.656),h=0.00698,phi=0.02536*(5.65/(0.656^2))),
  ## start = list(c=4,d=-0.03,phi=2,h=0.013),
  ## method="SANN",control=list(maxit=10000),
  data=bugdat)
g4_hyperq_bad = mle2(killed~dbinom(prob=1/(1/(c/(1+(size/d)+phi*(size/d)^2))+
                                 h*initial),
  size=initial),
  start=c(as.list(coef(g4_hyper)),list(phi=0)),
  data=bugdat)
g4_exp = mle2(killed~dbinom(prob=1/(1/(c*exp(-size/d))+
                                 h*initial),
  size=initial),
  start=sv,
  data=bugdat)
pred.p = predict(g4_powricker)/bugdat$initial
resids = with(bugdat,
  (predict(g4_powricker)-killed)/sqrt(initial*pred.p*(1-pred.p)))
xr = data.frame(bugdat,resid=resids)
disp = sum(resids^2)/(nrow(bugdat)-length(coef(g4_powricker)))
(q1 = ICtab(g4_ricker,g4_genricker,g4_powricker,g4_logist,g4_hyper,g4_hyperq,g4_exp,
   g4_hyperq_bad,
        type="qAIC",sort=TRUE,delta=TRUE,dispersion=disp,weights=TRUE))
bugmodlist = list(powricker=g4_powricker,
  genricker=g4_genricker,
  logist=g4_logist,
  hyperq=g4_hyperq,
  ricker=g4_ricker,
  exp=g4_exp,
  hyper=g4_hyper)
bugmodlist <- bugmodlist[gsub("g4_","",attr(q1,"row.names"))]
bugmodlist <- bugmodlist[!is.na(names(bugmodlist))]


###################################################
### chunk number 14: 
###################################################
#line 317 "mccoy-shape.Rnw"
predframe = expand.grid(size=seq(5,60,by=0.5),
  initial=sort(unique(x$initial)))
pfun = function(fit,name) {
  data.frame(predframe,killed=predict(fit,newdata=predframe),model=factor(name))
}
mframe= do.call(rbind,mapply(pfun,bugmodlist,names(bugmodlist),SIMPLIFY=FALSE))
belo_mframe = data.frame(mframe,predtype="belo")


###################################################
### chunk number 15: plot_models2
###################################################
#line 327 "mccoy-shape.Rnw"
print(ggplot(subset(belo_mframe,model!="hyper"),
             aes(x=size,y=killed,colour=factor(initial),
                 group=factor(initial)))+
      geom_line()+facet_wrap(~model)+geom_point(data=x)+
      labs(colour="initial density"))


###################################################
### chunk number 16: 
###################################################
#line 346 "mccoy-shape.Rnw"
ccvec = sapply(list(g3_genricker,g4_genricker,g3_powricker,g4_powricker),coef)
round(matrix(ccvec["gamma",],
       nrow=2,dimnames=list(c("dragonfly","waterbug"),
                c("genR","powR"))),3)


###################################################
### chunk number 17: 
###################################################
#line 363 "mccoy-shape.Rnw"
grc = coef(g3_genricker)
svec2 = c(as.list(grc[1:3]),
          list(gamma1=unname(grc["gamma"]),gamma2=unname(grc["gamma"])))
g3_genricker2 = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d)^gamma1*(exp(1-size/d))^gamma2)+
                       h*initial),
                     size=initial),
       start=svec2,
       method="Nelder-Mead",
  data=ododat)
grc = coef(g4_genricker)
svec2 = c(as.list(grc[1:3]),
          list(gamma1=unname(grc["gamma"]),gamma2=unname(grc["gamma"])))
g4_genricker2 = 
  mle2(killed~dbinom(prob=1/(1/(c*(size/d)^gamma1*(exp(1-size/d))^gamma2)+
                       h*initial),
                     size=initial),
       start=svec2,
       method="Nelder-Mead",
  data=bugdat)
ccvec = sapply(list(g3_genricker2,g4_genricker2),coef)
ccvec = t(ccvec[c("gamma1","gamma2"),])
rownames(ccvec) = c("dragonfly","waterbug")
round(ccvec,3)


###################################################
### chunk number 18: 
###################################################
#line 391 "mccoy-shape.Rnw"
h = g4_genricker2@details$hessian
round(eigen(h)$values,4)


###################################################
### chunk number 19: 
###################################################
#line 401 "mccoy-shape.Rnw"
mm = cbind(c(coef(g4_powricker),1),
  c(coef(g4_genricker),coef(g4_genricker)["gamma"]),coef(g4_genricker2))
mm = rbind(mm,
  c(-logLik(g4_powricker),-logLik(g4_genricker),-logLik(g4_genricker2)))
colnames(mm)=c("pow","gen","gen2")
rownames(mm)[4:6] = c("gamma1","gamma2","nlogLik")
round(mm,3)


###################################################
### chunk number 20: 
###################################################
#line 412 "mccoy-shape.Rnw"
mm = cbind(c(coef(g3_powricker),1),
  c(coef(g3_genricker),coef(g3_genricker)["gamma"]),coef(g3_genricker2))
colnames(mm)=c("pow","gen","gen2")
mm = rbind(mm,
  c(-logLik(g3_powricker),-logLik(g3_genricker),-logLik(g3_genricker2)))
rownames(mm)[4:6] = c("gamma1","gamma2","nlogLik")
round(mm,3)


###################################################
### chunk number 21:  eval=FALSE
###################################################
## #line 446 "mccoy-shape.Rnw"
## ## still working on it
## gvec = seq(1,4,length=20)
## res <- matrix(ncol=6,nrow=length(gvec)^2) 
## k <- 1  
## for (i in seq_along(gvec))
##   for (j in seq_along(gvec)) {
##     mtmp = try(mle2(killed~dbinom(prob=1/(1/(c*
##                                 (size/d)^gamma1*(exp(1-size/d))^gamma2)+
##                                 h*initial),
##       size=initial),
##       start=svec2,
##       method="Nelder-Mead",
##       data=bugdat,
##       fixed=list(gamma1=gvec[i],gamma2=gvec[j])))
##     if (inherits(mtmp,"try-error")) {
##       res[k,] <- rep(NA,6)
##     } else {
##       res[k,] <- c(-logLik(mtmp),coef(mtmp))
##     }
##     k <- k+1
## }  
## mm = matrix(res[,1],ncol=20)
## ## plot(row(mm),col(mm),type="n")
## ## text(row(mm),col(mm),round(mm-min(mm,na.rm=TRUE)+1,2),cex=0.5)
## apply(mm,2,sd,na.rm=TRUE)


###################################################
### chunk number 22: 
###################################################
#line 474 "mccoy-shape.Rnw"
all_mframe=rbind(odo_mframe,belo_mframe)
mframe_powricker=subset(all_mframe,model="powricker")


###################################################
### chunk number 23:  eval=FALSE
###################################################
## #line 488 "mccoy-shape.Rnw"
## download.file("http://admb-project.googlecode.com/files/admb-ide-345-6.exe")
## system("admb-ide-345-6.exe")


###################################################
### chunk number 24: 
###################################################
#line 493 "mccoy-shape.Rnw"
source("admb-funs.R")
admb_setup()
## input data
inp =   c(list(nobs = nrow(ododat)), 
  subset(ododat, select = c(killed,
                   size, initial)))



###################################################
### chunk number 25: 
###################################################
#line 503 "mccoy-shape.Rnw"
d3 = do_admb("mccoypred2",inp,c(sv,list(g=1)))
round(cbind(coef(d3),coef(g3_powricker)),4)
c(logLik(d3)-logLik(g3_powricker))


###################################################
### chunk number 26: admbreff
###################################################
#line 522 "mccoy-shape.Rnw"
nblock = length(levels(ododat$block))
inp2 = c(list(nobs=nrow(ododat),
  nblock=nblock),
  subset(ododat,select=c(killed,size,initial)),
  list(Z=model.matrix(~block-1,ododat)))
pars = c(as.list(coef(d3)),list(sigma_c=0.05,
  u=rep(0,nblock)))
t1.admb <- system.time(d4 <- do_admb("mccoypred6",inp2,pars,re=TRUE))
fn <- "d4_admb.RData"
if (!file.exists(fn)) {
  t2.admb <- system.time(d4m <- do_admb("mccoypred6",inp2,pars,re=TRUE,
                                        safe=FALSE,
                                        mcmc=TRUE))
  t3.admb <- system.time(d4m2 <- do_admb("mccoypred6",inp2,pars,re=TRUE,
                                         safe=FALSE,
                                         mcmc=TRUE,mcmcsteps=10000))
  save("t2.admb","d4m","t3.admb","d4m2",file=fn)
} else load(fn)


###################################################
### chunk number 27: 
###################################################
#line 544 "mccoy-shape.Rnw"
mm = cbind(coef(d4)[1:4],coef(g3_powricker))
colnames(mm)=c("ADMB","base")
round(mm,3)


###################################################
### chunk number 28: 
###################################################
#line 551 "mccoy-shape.Rnw"
with(as.list(coef(d4)),sigma_c/c)


###################################################
### chunk number 29: 
###################################################
#line 559 "mccoy-shape.Rnw"
logLik(d4)-logLik(d3)


###################################################
### chunk number 30: 
###################################################
#line 569 "mccoy-shape.Rnw"
plot(d4m$hist,pars=1:5)
print(last_plot())


###################################################
### chunk number 31: 
###################################################
#line 581 "mccoy-shape.Rnw"
library(R2WinBUGS)
xstart = c(as.list(coef(g3_powricker)),list(sd.c=0.1))
set.seed(1001)
xstart2 = list(
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)),
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)),
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)))
xx <- ododat
xx$block = as.numeric(xx$block)
dat = c(as.list(subset(xx,select=c(killed,size,initial,block))),
  N=nrow(ododat),
  nblock=length(levels(ododat$block)))


###################################################
### chunk number 32: 
###################################################
#line 603 "mccoy-shape.Rnw"
fn = "mccoybugs1.RData"
if (!file.exists(fn)) {
  t1.wbugs <- system.time(b1 <- bugs(data=dat,
    inits=c(list(xstart),xstart2),
    n.chains=4,
    n.iter=80000,
    parameters=c(names(xstart),"cvec"),
    model.file="mccoypred.bugs",
    working.directory=getwd()))
  save("b1","t1.wbugs",file=fn)
  } else load(fn)
fn <- "" ## prevent overwriting later on


###################################################
### chunk number 33: 
###################################################
#line 618 "mccoy-shape.Rnw"
g3prof_powricker = profile(g3_powricker_pen)


###################################################
### chunk number 34: 
###################################################
#line 637 "mccoy-shape.Rnw"
pp = plot(d4m2$hist,pars=1:5)
ppp = read_psv("mccoypred6")
mpp = melt(ppp[1:5])
names(mpp)[1] = "param"
zz = ddply(mpp,"param",function(x) {dd <- density(x$value);
                                    data.frame(X1=dd$x,X2=dd$y)})
pp+geom_line(data=zz,colour="red")
library(emdbook)
b2 = lump.mcmc.list(as.mcmc(b1))
mpp2 = melt(as.data.frame(b2)[1:5])
names(mpp2)[1] = "param"
levels(mpp2$param) <- c("rc","rd","rh","rg","rsigma_c")
zz2 = ddply(mpp2,"param",function(x) {dd <- density(x$value);
                                    data.frame(X1=dd$x,X2=dd$y)})
print(pp+geom_line(data=zz,colour="red")+
  geom_line(data=zz2,colour="blue"))


###################################################
### chunk number 35: 
###################################################
#line 660 "mccoy-shape.Rnw"
ci1 = melt(HPDinterval(b2[,1:5]))
ci2 = melt(HPDinterval(as.mcmc(ppp[1:5])))
levels(ci1$X1) <- levels(ci2$X1)
ci3 = melt(confint(g3prof_powricker))
levels(ci3$X2) <- levels(ci2$X2)
levels(ci3$X1) <- levels(ci2$X1)[1:4]
ciframe = data.frame(do.call(rbind,list(ci1,ci2,ci3)),
  method=factor(rep(c("WinBUGS","ADMB","mle2"),c(10,10,8))))
est_mle = data.frame(method="mle2",X2="est",
  value=coef(g3_powricker),X1=names(coef(g3_powricker)))
levels(est_mle$X1)=levels(ci3$X1)
est_ADMB = data.frame(method="ADMB",X2="est",
  value=coef(d4m2)[1:5],
  X1=names(coef(d4m2))[1:5])
levels(est_ADMB$X1)=levels(ci2$X1)
est_WinBUGS = data.frame(method="WinBUGS",X2="est",
  value=unlist(b1$mean[1:5]),
  X1=names(b1$mean)[1:5])
levels(est_WinBUGS$X1)=levels(ci2$X1)
estframe = rbind(est_ADMB,est_mle,est_WinBUGS)
ciframe = rbind(ciframe,estframe)
ciframe = recast(ciframe,...~X2)
print(ggplot(ciframe,aes(x=method,y=est,ymin=lower,ymax=upper))+geom_point()+
       geom_errorbar()+facet_wrap(~X1,scale="free"))


###################################################
### chunk number 36: 
###################################################
#line 689 "mccoy-shape.Rnw"
dd2 = density(b1$sims.list$sd.c)
matrix(c(b1$mean$sd.c,
         dd2$x[which.max(dd2$y)],
         mean(ppp$rsigma_c),
         coef(d4m2)["sigma_c"]),
       nrow=2,
       dimnames=list(c("mean","mode"),
         c("WinBUGS","ADMB")))                        


###################################################
### chunk number 37:  eval=FALSE
###################################################
## #line 701 "mccoy-shape.Rnw"
## xyplot(as.mcmc(ppp[1:5]),layout=NULL)
## xyplot(as.mcmc(b1),layout=NULL)


