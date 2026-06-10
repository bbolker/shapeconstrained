x = read.csv("Response surfaces data_Gamboa _data.csv")
library(gdata)
x$block = factor(x$block)
x <- drop.levels(subset(x,cohort=="single"))

load("modelfits.RData")

library(ggplot2)
ss = scale_shape_manual(name="block",values=1:7)
ss2A =   scale_linetype_manual("size class",
                        c("21","31","42","53","74"))
ss2B =   scale_linetype_manual("initial density",
  c("21","31","42","53","63","74","84"))

lfun = function(variable,value) {
  ifelse(value=="belo","Water bug","Dragonfly")
}

## figure out size classes, again
with(x,tapply(size,list(sizeclass),mean))
## 8, 9.5, 18, 30, 45?
library(gdata)
mframe = drop.levels(subset(all_mframe,model=="powricker" &
  size %in% c(8,9.5,18,30,45)))
mframe$sizeclass = factor(mframe$size)
levels(mframe$sizeclass) <- as.character(as.numeric(mframe$sizeclass))

## take out block info -- too messy
p1 = ggplot(x,aes(x=initial,y=killed,
  group=factor(sizeclass),
  lty=factor(sizeclass),
  pch=factor(sizeclass),
  fill=factor(sizeclass)))
## currently unused:   ## colour=factor(sizeclass)))+
p2 <- p1 +
  scale_fill_grey(start=1,end=0)+
  scale_shape_manual("size class",value=21:25)+
  ## hack: avoid line type 3, problems with PDF output
  ss2A+
  labs(fill="size class",
       lty="size class")+
  facet_grid(~predtype,labeller=lfun)+ ## split by pred type
  stat_sum(aes(size=..n..))+  ## point size according to overlap
  scale_size_continuous(to=c(2,6),legend=FALSE)+ ## scrunch size range of points
  theme_bw(base_size=16)+ ## black & white; increase font size etc.
  xlab("Initial density (#/tank)")+ylab("Number killed") ## axis labels

p3 = p2 + geom_line(data=mframe)


mframe = drop.levels(subset(all_mframe,model=="powricker"))

p4 = ggplot(x,aes(x=size,y=killed,
  group=factor(initial),
  lty=factor(initial),
  pch=factor(initial),
  fill=factor(initial)))
## currently unused:   ## colour=factor(sizeclass)))+
p5 <- p4 +
  scale_fill_grey(start=1,end=0)+
  scale_shape_manual("initial density",value=c(21:25,7:8))+
  ss2B+
  labs(fill="initial density",lty="initial density")+
  facet_grid(~predtype,labeller=lfun)+ ## split by pred type
  geom_point() +  ## point size according to overlap
### scale_size_continuous(to=c(2,6),legend=FALSE)+ ## scrunch size range of points
  theme_bw(base_size=16)+ ## black & white; increase font size etc.
  xlab("Size (mm)")+ylab("Number killed") ## axis labels

p6 = p5 + geom_line(data=mframe)

vptop <- viewport(width = 0.95, height = 0.5, x = 0.475, y = 0.75)
vpbottom <- viewport(width = 1.0, height = 0.5, x=0.5, y=0.25)

pdf("mwm_modelfit1.pdf")
## ,height=5,width=10)
print(p3, vp=vptop)
print(p6, vp=vpbottom)
##  geom_line(data=mframe))
dev.off()



