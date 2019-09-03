#clean up working environment, define plotting parameters, set default palette...
rm(list=ls())
def.par<-par(no.readonly=T)
alter.cols<-function(x,alph=0.5,darken=1,name=F){
  if(name){
    cols<-col2rgb(x,alpha=T)/255
  }else{
    cols<-col2rgb(as.numeric(sort(x)),alpha=T)/255
  }
  cols[1:3,]<-cols[1:3,]/darken
  cols[4,]<-alph
  return(mapply(rgb,red=cols[1,],green=cols[2,],blue=cols[3,],alpha=cols[4,]))
}
custom.palette<-c("red","orange","green2","blue","purple")
palette(colorRampPalette(custom.palette)(10))
library(abind)

#Load up auxillary functions
source("aux_fun/GeomFunctions.R")
source("aux_fun/DataWrangFunctions.R")

#load data, replace empty cells w/ NA's
raw<-read.csv("data/vib_lms.csv",na.strings="")
#trim data
all.nas<-function(x){
  all(is.na(x))
}
raw<-raw[,-which(apply(raw,2,all.nas))]

#generate initial covariance matrix and landmark array
covars<-gen.covar.mat(raw,covar.start=3,lm.start=5)
lms<-arrange.data(raw,lm.start=5)

#split leaves into left and right halves
splits<-split.leaves(lms,covars,EFN.info.start=4)

#transform landmarks using procrustes procedure
trans.lms<-proc.fit(splits$x,scale=T)
trans.lms[grep("EFN14",dimnames(trans.lms)[[1]]),,"168.R"]<-NA

#plotting whole leaf landmarks to make sure code worked!
#plot(trans.lms[1:3,"y",]~trans.lms[1:3,"x",],col=rep(c("red","blue","green"),dim(trans.lms)[3]),asp=1)
#segments(x0=trans.lms[1:3,"x",],y0=trans.lms[1:3,"y",],
         #x1=trans.lms[c(2,3,1),"x",],y1=trans.lms[c(2,3,1),"y",],col=rgb(0,0,0,0.1))
#text(x=trans.lms[1:3,"x",],y=trans.lms[1:3,"y",],labels=rep(dimnames(trans.lms)[[3]],each=3))

#EFN locations
EFNs.max.1<-extract.points(trans.lms,"EFN\\d+.max.1")
EFNs.max.2<-extract.points(trans.lms,"EFN\\d+.max.2")
EFNs.min.1<-extract.points(trans.lms,"EFN\\d+.min.1")
EFNs.min.2<-extract.points(trans.lms,"EFN\\d+.min.2")
EFNs.base<-extract.points(trans.lms,"EFN\\d+.base")
EFNs.tip<-extract.points(trans.lms,"EFN\\d+.tip")
EFNs.ys<-colMeans(rbind(extract.ys(EFNs.max.1),extract.ys(EFNs.max.2)))
EFNs.xs<-colMeans(rbind(extract.xs(EFNs.max.1),extract.xs(EFNs.max.2)))
EFNs.ys2<-colMeans(rbind(extract.ys(EFNs.base),extract.ys(EFNs.tip)))
EFNs.xs2<-colMeans(rbind(extract.xs(EFNs.base),extract.xs(EFNs.tip)))
side.view.EFNs<-which(is.na(EFNs.xs)&!is.na(EFNs.xs2));names(side.view.EFNs)<-NULL
EFNs.xs[side.view.EFNs]<-EFNs.xs2[side.view.EFNs];EFNs.ys[side.view.EFNs]<-EFNs.ys2[side.view.EFNs]
d2marg<-mapply(dist.p,p1=split(cbind(EFNs.xs,EFNs.ys),1:length(EFNs.xs)),
               p2=extract.points(trans.lms,"EFN\\d+.margin"))
d2vein<-mapply(dist.p,p1=split(cbind(EFNs.xs,EFNs.ys),1:length(EFNs.xs)),
               p2=extract.points(trans.lms,"EFN\\d+.vein"))
vein2marg<-mapply(dist.p,p1=extract.points(trans.lms,"EFN\\d+.vein"),
                  p2=extract.points(trans.lms,"EFN\\d+.marg"))
d2marg<-ifelse(d2vein>vein2marg,-1*d2marg,d2marg)
names(d2marg)<-names(EFNs.xs)

trans.lms<-proc.fit(splits$x,scale=F)
trans.lms[grep("EFN14",dimnames(trans.lms)[[1]]),,"168.R"]<-NA
#EFN areas/widths
EFNs.max.widths<-mapply(dist.p,p1=extract.points(trans.lms,"EFN\\d+.max.1"),
                       p2=extract.points(trans.lms,"EFN\\d+.max.2"))
EFNs.min.widths<-mapply(dist.p,p1=extract.points(trans.lms,"EFN\\d+.min.1"),
                       p2=extract.points(trans.lms,"EFN\\d+.min.2"))
EFNs.max.widths.new<-ifelse(EFNs.max.widths>EFNs.min.widths,EFNs.max.widths,EFNs.min.widths)
EFNs.min.widths.new<-ifelse(EFNs.min.widths<EFNs.max.widths,EFNs.min.widths,EFNs.max.widths)
#the test will reutrn na with missing values; this results in loss of info, so exclude returned na's from replacement
EFNs.max.widths[!is.na(EFNs.max.widths.new)]<-EFNs.max.widths.new[!is.na(EFNs.max.widths.new)]
EFNs.min.widths[!is.na(EFNs.min.widths.new)]<-EFNs.min.widths.new[!is.na(EFNs.min.widths.new)]
EFNs.areas<-pi*EFNs.max.widths/2*EFNs.min.widths/2

EFNs.max.bowl.widths<-mapply(dist.p,p1=extract.points(trans.lms,"EFN\\d+.bowl.max.1"),
                            p2=extract.points(trans.lms,"EFN\\d+.bowl.max.2"))
EFNs.min.bowl.widths<-mapply(dist.p,p1=extract.points(trans.lms,"EFN\\d+.bowl.min.1"),
                            p2=extract.points(trans.lms,"EFN\\d+.bowl.min.2"))
EFNs.max.bowl.widths.new<-ifelse(EFNs.max.bowl.widths>EFNs.min.bowl.widths,EFNs.max.bowl.widths,EFNs.min.bowl.widths)
EFNs.min.bowl.widths.new<-ifelse(EFNs.min.bowl.widths<EFNs.max.bowl.widths,EFNs.min.bowl.widths,EFNs.max.bowl.widths)
EFNs.max.bowl.widths[!is.na(EFNs.max.bowl.widths.new)]<-EFNs.max.bowl.widths.new[!is.na(EFNs.max.bowl.widths.new)]
EFNs.min.bowl.widths[!is.na(EFNs.min.bowl.widths.new)]<-EFNs.min.bowl.widths.new[!is.na(EFNs.min.bowl.widths.new)]
EFNs.bowl.areas<-pi*EFNs.max.bowl.widths/2*EFNs.min.bowl.widths/2

EFNs.heights<-mapply(dist.p,p1=extract.points(trans.lms,"EFN\\d+.base"),
                     p2=extract.points(trans.lms,"EFN\\d+.tip"))

#getting species factor vector
sp_fac<-spec.var(EFNs.xs,splits$x_covar,"sp")

spec_fac<-sapply(strsplit(names(EFNs.areas),' '),'[',1)

spec_length<-mapply(dist.p,p1=extract.points(trans.lms,"^tip$"),
                    p2=extract.points(trans.lms,"^base$"))
names(spec_length)<-sapply(strsplit(names(spec_length),' '),'[',1)
spec_length<-spec_length[spec_fac]

#getting sp-specific mean, std, and se measurements (right now for EFN areas)
get.sp.means<-function(trait,cov.mat,log=T){
  spec_fac<-sapply(strsplit(names(trait),'\\.'),'[',1)
  if(log){
    spec.means<-tapply(log(trait),spec_fac,mean,na.rm=T)
  }else{
    spec.means<-tapply(trait,spec_fac,mean,na.rm=T)
  }
  tmp.sp_fac<-spec.var(spec.means,cov.mat,"sp")
  sp.means<-tapply(spec.means,tmp.sp_fac,mean,na.rm=T)
  sp.sd<-tapply(spec.means,tmp.sp_fac,sd,na.rm=T)
  sp.se<-tapply(spec.means,tmp.sp_fac,function(x,...) {sd(x,...)/sqrt(length(!is.na(x)))},na.rm=T)
  return(list(mean=sp.means,sd=sp.sd,se=sp.se))
}

plot.sim<-function(trait,sp.means,cov.mat,y.label,
                   sd=T,sim.log=T,n.sim=20,jitter.fac=1,trans.fac=0.15,sim.trans.fac=0.2,sim.cex=2,name.cex=0.6){
  if(sd){
    error<-sp.means$sd
  }else{
    error<-sp.means$se
  }
  tmp.sp_fac<-spec.var(trait,cov.mat,"sp")
  plot(trait[order(tmp.sp_fac)]~jitter(as.numeric(sort(tmp.sp_fac)),jitter.fac),
       col=alter.cols(tmp.sp_fac,trans.fac),pch=1,xaxt='n',xlab='',ylab=y.label)
  if(sim.log){
    points(exp(rnorm(length(sp.means$mean)*n.sim,rep(sp.means$mean,each=n.sim),rep(error,each=n.sim)))~
             rep(1:length(sp.means$mean),each=n.sim),
           col=alter.cols(rep(1:length(sp.means$mean),each=n.sim),sim.trans.fac),pch=16,cex=sim.cex)
  }else{
    points(rnorm(length(sp.means$mean)*n.sim,rep(sp.means$mean,each=n.sim),rep(error,each=n.sim))~
             rep(1:length(sp.means$mean),each=n.sim),
           col=alter.cols(rep(1:length(sp.means$mean),each=n.sim),sim.trans.fac),pch=16,cex=sim.cex)
  }
  axis(side=1,labels=sort(unique(tmp.sp_fac)),at=1:length(unique(tmp.sp_fac)),las=2,cex.axis=0.6)
}

sp.areas<-get.sp.means(EFNs.areas,covars)
plot.sim(EFNs.areas,sp.areas,splits$x_covar,"EFN area (mm)")


trans.lms<-proc.fit(splits$x,scale=T)
trans.lms[grep("EFN14",dimnames(trans.lms)[[1]]),,"168.R"]<-NA
#Fitting beta distribution to x and y distributions of EFNs across unit leaf scale...
library(fitdistrplus)
std.EFNs.xs<-EFNs.xs/trans.lms["side.max","x",sapply(strsplit(names(EFNs.areas),' '),'[',1)]
std.EFNs.xs[std.EFNs.xs< -0.1]<-NA
midrib.pos<-0-min(std.EFNs.xs,na.rm=T)
marg.pos<-1-min(std.EFNs.xs,na.rm=T)
std.EFNs.xs<-std.EFNs.xs-min(std.EFNs.xs,na.rm=T)
midrib.pos<-midrib.pos/max(marg.pos,std.EFNs.xs,na.rm=T)
marg.pos<-marg.pos/max(marg.pos,std.EFNs.xs,na.rm=T)
std.EFNs.xs<-std.EFNs.xs/max(marg.pos,std.EFNs.xs,na.rm=T)
std.EFNs.ys<-EFNs.ys/trans.lms["tip","y",sapply(strsplit(names(EFNs.areas),' '),'[',1)]
base.pos<-0-min(std.EFNs.ys,na.rm=T)
tip.pos<-1-min(std.EFNs.ys,na.rm=T)
std.EFNs.ys<-std.EFNs.ys-min(std.EFNs.ys,na.rm=T)
base.pos<-base.pos/max(tip.pos,std.EFNs.ys,na.rm=T)
tip.pos<-tip.pos/max(tip.pos,std.EFNs.ys,na.rm=T)
std.EFNs.ys<-std.EFNs.ys/max(tip.pos,std.EFNs.ys,na.rm=T)

sp.EFNs.xs<-tapply(std.EFNs.xs[!is.na(std.EFNs.xs)],sp_fac[!is.na(std.EFNs.xs)],fitdist,dist="beta",method="mme")
sp.EFNs.ys<-tapply(std.EFNs.ys[!is.na(std.EFNs.ys)],sp_fac[!is.na(std.EFNs.ys)],fitdist,dist="beta",method="mme")

sp.EFNs.xs.shape1<-sapply(sapply(sp.EFNs.xs,'[',1),'[',1)
names(sp.EFNs.xs.shape1)<-sapply(strsplit(names(sp.EFNs.xs.shape1),"\\."),'[',1)
sp.EFNs.xs.shape2<-sapply(sapply(sp.EFNs.xs,'[',1),'[',2)
names(sp.EFNs.xs.shape2)<-sapply(strsplit(names(sp.EFNs.xs.shape2),"\\."),'[',1)
sp.EFNs.ys.shape1<-sapply(sapply(sp.EFNs.ys,'[',1),'[',1)
names(sp.EFNs.ys.shape1)<-sapply(strsplit(names(sp.EFNs.ys.shape1),"\\."),'[',1)
sp.EFNs.ys.shape2<-sapply(sapply(sp.EFNs.ys,'[',1),'[',2)
names(sp.EFNs.ys.shape2)<-sapply(strsplit(names(sp.EFNs.ys.shape2),"\\."),'[',1)

EFNs.xs.mus<-sp.EFNs.xs.shape1/(sp.EFNs.xs.shape1+sp.EFNs.xs.shape2)
EFNs.xs.rates<-1/(sp.EFNs.xs.shape1+sp.EFNs.xs.shape2)

EFNs.ys.mus<-sp.EFNs.ys.shape1/(sp.EFNs.ys.shape1+sp.EFNs.ys.shape2)
EFNs.ys.rates<-1/(sp.EFNs.ys.shape1+sp.EFNs.ys.shape2)


#bootstrap
iter<-100
perm.EFNs.xs.mus<-matrix(rep(NA),nrow=length(unique(sp_fac)),ncol=iter)
rownames(perm.EFNs.xs.mus)<-sort(as.character(unique(sp_fac)))
perm.EFNs.xs.rates<-matrix(rep(NA),nrow=length(unique(sp_fac)),ncol=iter)
rownames(perm.EFNs.xs.rates)<-sort(as.character(unique(sp_fac)))
perm.EFNs.ys.mus<-matrix(rep(NA),nrow=length(unique(sp_fac)),ncol=iter)
rownames(perm.EFNs.ys.mus)<-sort(as.character(unique(sp_fac)))
perm.EFNs.ys.rates<-matrix(rep(NA),nrow=length(unique(sp_fac)),ncol=iter)
rownames(perm.EFNs.ys.rates)<-sort(as.character(unique(sp_fac)))
perm.EFNs.xs<-std.EFNs.xs
perm.EFNs.ys<-std.EFNs.ys
for(i in 1:iter){
  for(j in sort(as.character(unique(sp_fac)))){
    perm.EFNs.xs[!(is.na(std.EFNs.xs))&sp_fac==j]<-sample(std.EFNs.xs[!(is.na(std.EFNs.xs))&sp_fac==j],replace=T)
    perm.EFNs.ys[!(is.na(std.EFNs.ys))&sp_fac==j]<-sample(std.EFNs.xs[!(is.na(std.EFNs.ys))&sp_fac==j],replace=T)
  }
  perm.sp.EFNs.xs<-tapply(perm.EFNs.xs[!is.na(perm.EFNs.xs)],
                          sp_fac[!is.na(perm.EFNs.xs)],fitdist,dist="beta",method="mme")
  perm.sp.EFNs.ys<-tapply(perm.EFNs.ys[!is.na(perm.EFNs.ys)],
                          sp_fac[!is.na(perm.EFNs.ys)],fitdist,dist="beta",method="mme")
  perm.sp.EFNs.xs.shape1<-sapply(sapply(perm.sp.EFNs.xs,'[',1),'[',1)
  names(perm.sp.EFNs.xs.shape1)<-sapply(strsplit(names(perm.sp.EFNs.xs.shape1),"\\."),'[',1)
  perm.sp.EFNs.xs.shape2<-sapply(sapply(perm.sp.EFNs.xs,'[',1),'[',2)
  names(perm.sp.EFNs.xs.shape2)<-sapply(strsplit(names(perm.sp.EFNs.xs.shape2),"\\."),'[',1)
  perm.sp.EFNs.ys.shape1<-sapply(sapply(perm.sp.EFNs.ys,'[',1),'[',1)
  names(perm.sp.EFNs.ys.shape1)<-sapply(strsplit(names(perm.sp.EFNs.ys.shape1),"\\."),'[',1)
  perm.sp.EFNs.ys.shape2<-sapply(sapply(perm.sp.EFNs.ys,'[',1),'[',2)
  names(perm.sp.EFNs.ys.shape2)<-sapply(strsplit(names(perm.sp.EFNs.ys.shape2),"\\."),'[',1)
  perm.EFNs.xs.mus[,i]<-perm.sp.EFNs.xs.shape1/(perm.sp.EFNs.xs.shape1+perm.sp.EFNs.xs.shape2)
  perm.EFNs.xs.rates[,i]<-1/(perm.sp.EFNs.xs.shape1+perm.sp.EFNs.xs.shape2)
  perm.EFNs.ys.mus[,i]<-perm.sp.EFNs.ys.shape1/(perm.sp.EFNs.ys.shape1+perm.sp.EFNs.ys.shape2)
  perm.EFNs.ys.rates[,i]<-1/(perm.sp.EFNs.ys.shape1+perm.sp.EFNs.ys.shape2)
}

EFNs.xs.mus.se<-EFNs.xs.mus
EFNs.xs.rates.se<-EFNs.xs.mus
EFNs.ys.mus.se<-EFNs.xs.mus
EFNs.ys.rates.se<-EFNs.xs.mus
for(i in sort(as.character(unique(sp_fac)))){
  EFNs.xs.mus.se[i]<-sd(perm.EFNs.xs.mus[i,])
  EFNs.xs.rates.se[i]<-sd(perm.EFNs.xs.rates[i,])
  EFNs.ys.mus.se[i]<-sd(perm.EFNs.ys.mus[i,])
  EFNs.ys.rates.se[i]<-sd(perm.EFNs.ys.rates[i,])
}

#get 95% errors (not assuming normal, not worrying se vs. sd distinction)
#EFNs.xs.mus.95es<-EFNs.xs.mus
#EFNs.xs.rates.95es<-EFNs.xs.mus
#EFNs.ys.mus.95es<-EFNs.xs.mus
#EFNs.ys.rates.95es<-EFNs.xs.mus
#for(i in sort(as.character(unique(sp_fac)))){
  #EFNs.xs.mus.95es[i]<-mean(abs(EFNs.xs.mus[i]-quantile(perm.EFNs.xs.mus[i,],c(0.025,0.975))))
  #EFNs.xs.rates.95es[i]<-mean(abs(EFNs.xs.rates[i]-quantile(perm.EFNs.xs.rates[i,],c(0.025,0.975))))
  #EFNs.ys.mus.95es[i]<-mean(abs(EFNs.ys.mus[i]-quantile(perm.EFNs.ys.mus[i,],c(0.025,0.975))))
  #EFNs.ys.rates.95es[i]<-mean(abs(EFNs.ys.rates[i]-quantile(perm.EFNs.ys.rates[i,],c(0.025,0.975))))
#}


plot(std.EFNs.ys[sp_fac=="hengshanicum"]~std.EFNs.xs[sp_fac=="hengshanicum"])

plot(EFNs.ys.mus)

#Fitting poisson params to EFN counts
EFNs.obs<-extract.xs(extract.points(trans.lms,"EFN\\d+.max.1"))
names(EFNs.obs)<-sapply(strsplit(sapply(strsplit(names(EFNs.obs)," "),'[',1),"\\."),'[',1)
EFNs.obs<-ifelse(is.na(EFNs.obs),F,T)
spec.EFNs.nos<-tapply(EFNs.obs,spec_fac,sum)
half.specs<-tapply(names(EFNs.obs),names(EFNs.obs),length)
half.specs<-names(half.specs)[half.specs==min(half.specs)]
spec.EFNs.nos[names(spec.EFNs.nos)%in%half.specs]<-2*spec.EFNs.nos[names(spec.EFNs.nos)%in%half.specs]

no.specs<-tapply(spec.EFNs.nos,sp_fac2,length)
one.spec<-names(no.specs)[no.specs==1]
sp.EFNs.nos.fit<-tapply(spec.EFNs.nos[!(sp_fac2%in%one.spec)],
                    sp_fac2[!(sp_fac2%in%one.spec)],fitdist,dist="pois",method="mle")
sp.EFNs.nos<-unlist(sapply(sp.EFNs.nos.fit,'[',1))
names(sp.EFNs.nos)<-sapply(strsplit(names(sp.EFNs.nos),"\\."),'[',1)
one.spec.nos<-spec.EFNs.nos[sp_fac2%in%one.spec]
one.spec.nos<-one.spec.nos[order(as.numeric(names(one.spec.nos)))]
names(one.spec.nos)<-one.spec
sp.EFNs.nos<-c(sp.EFNs.nos,one.spec.nos)[sort(as.character(unique(sp_fac)))]

sp.EFNs.nos.se<-unlist(sapply(sp.EFNs.nos.fit,'[',3))
names(sp.EFNs.nos.se)<-sapply(strsplit(names(sp.EFNs.nos.se),"\\."),'[',1)
one.spec.nos.se<-rep(NA,length(one.spec))
names(one.spec.nos.se)<-one.spec
sp.EFNs.nos.se<-c(sp.EFNs.nos.se,one.spec.nos.se)[sort(as.character(unique(sp_fac)))]

plot(spec.EFNs.nos[order(sp_fac2)]~jitter(as.numeric(sort(sp_fac2)),1),col=alter.cols(sp_fac2,0.25),pch=1,
     xaxt='n',xlab='',ylab="Number of EFNs")
n.sim<-50
points(rpois(length(sp.EFNs.nos)*n.sim,rep(sp.EFNs.nos,each=n.sim)+
               rnorm(length(sp.EFNs.nos)*n.sim,rep(sp.EFNs.nos.se,each=n.sim)))~
         jitter(rep(1:length(sp.EFNs.nos),each=n.sim),1),
       col=alter.cols(rep(1:length(sp.EFNs.nos),each=n.sim),0.05),pch=16,cex=2)
axis(side=1,labels=sort(unique(sp_fac)),at=1:length(unique(sp_fac)),las=2,cex.axis=0.6)



sp.min.widths<-get.sp.means(EFNs.min.widths,covars)
plot.sim(EFNs.min.widths,sp.min.widths,splits$x_covar,"EFN minimum width (mm)")

sp.max.widths<-get.sp.means(EFNs.max.widths,covars)
plot.sim(EFNs.max.widths,sp.max.widths,splits$x_covar,"EFN maximum width (mm)")

sp.min.bowl.widths<-get.sp.means(EFNs.min.bowl.widths,covars)
plot.sim(EFNs.min.bowl.widths,sp.min.bowl.widths,splits$x_covar,"EFN minimum secretory surface width (mm)")

sp.max.bowl.widths<-get.sp.means(EFNs.max.bowl.widths,covars)
plot.sim(EFNs.max.bowl.widths,sp.max.bowl.widths,splits$x_covar,"EFN maximum secretory surface width (mm)")

sp.heights<-get.sp.means(EFNs.heights,covars)
plot.sim(EFNs.heights,sp.heights,splits$x_covar,"EFN height (mm)")

marg.pos2<- -min(d2marg,na.rm=T)
std.d2marg<-d2marg-min(d2marg,na.rm=T)
sp.d2marg<-get.sp.means(std.d2marg,covars)
plot.sim(std.d2marg,sp.d2marg,splits$x_covar,"distance to margin")
abline(h=marg.pos2)

library(plotrix)
x.adv<-1.1
y.adv<-0.4
n.spp<-length(unique(sp_fac))
plot(0,xlim=c(-x.adv,n.spp*x.adv),ylim=c(-y.adv,n.spp*y.adv),col="white",asp=1)
draw.ellipse(x=seq(0,(n.spp-1)*x.adv,x.adv),y=seq(0,(n.spp-1)*y.adv,y.adv),
             a=exp(sp.max.widths$mean)/2,b=exp(sp.min.widths$mean)/2,angle=rep(0,n.spp),
             col="chartreuse3",border=NA)
draw.ellipse(x=seq(0,(n.spp-1)*x.adv,x.adv),y=seq(0,(n.spp-1)*y.adv,y.adv),
             a=exp(sp.max.bowl.widths$mean)/2,b=exp(sp.min.bowl.widths$mean)/2,angle=rep(0,n.spp),
             col="darkorange3",border=NA)
text(x=seq(0,(n.spp-1)*x.adv,x.adv)+rep(c(1,-1))*exp(sp.max.widths$mean)/2,
     y=seq(0,(n.spp-1)*y.adv,y.adv)+rep(c(-1,1))*exp(sp.min.widths$mean)/2,pos=rep(c(1,3)),
     names(sp.max.widths$mean),cex=0.6)

plot(0,xlim=c(-x.adv,n.spp*x.adv),ylim=c(-y.adv,n.spp*y.adv),col="white",asp=1)
for(i in 1:n.spp){
  xx<-(i-1)*x.adv
  yy<-(i-1)*y.adv
  polygon(x=c(xx-exp(sp.max.widths$mean)[i]/2,xx-exp(sp.max.bowl.widths$mean)[i]/2,
              xx+exp(sp.max.bowl.widths$mean)[i]/2,xx+exp(sp.max.widths$mean)[i]/2),
          y=c(yy,yy+exp(sp.heights$mean)[i],yy+exp(sp.heights$mean)[i],yy),col="chartreuse3",border=NA)
  draw.ellipse(x=xx,y=yy+exp(sp.heights$mean)[i],
               a=exp(sp.max.bowl.widths$mean)[i]/2,b=0.2*exp(sp.heights$mean)[i],angle=0,
               col="darkorange3",border=NA)
}
text(x=seq(0,(n.spp-1)*x.adv,x.adv),
     y=seq(0,(n.spp-1)*y.adv,y.adv)+rep(c(0,1))*exp(sp.heights$mean),pos=rep(c(1,3)),
     names(sp.max.widths$mean),cex=0.6)

###climate stuff###
climate.means<-read.csv("data/Viburnum.climateMeans.csv",row.names=1)
rownames(climate.means)<-sapply(strsplit(as.character(climate.means[,1])," "),'[',2)
climate.means<-climate.means[,-c(1:2)]
plot(climate.means[,"bio4"])

climate.pruned<-climate.means[rownames(climate.means)%in%pruned.tree$tip.label,]
plot(climate.pruned[,"bio4"]~climate.pruned[,"lat"])
text(x=climate.pruned[,"lat"],y=climate.pruned[,"bio4"],labels=rownames(climate.pruned))

sp.lats<-abs(climate.pruned[,"lat"]);names(sp.lats)<-rownames(climate.pruned)
sp.seasonality<-climate.pruned[,"bio4"];names(sp.seasonality)<-rownames(climate.pruned)
sp.temp2<-ifelse(sp.lats>23.5,0,1)

#exporting trait data
trait.dat<-cbind(EFNs.xs.mus,EFNs.xs.rates,EFNs.ys.mus,EFNs.ys.rates,
                 sp.max.widths$mean,sp.min.widths$mean,sp.max.bowl.widths$mean,sp.min.bowl.widths$mean,sp.heights$mean)
colnames(trait.dat)<-c("x.cent","xs.disp","ys.cent","ys.disp",
                       "max.width","min.width","ss.max.width","ss.min.width","height")
trait.dat[is.nan(trait.dat)]<-NA
trait.se<-cbind(EFNs.xs.mus.se,EFNs.xs.rates.se,EFNs.ys.mus.se,EFNs.ys.rates.se,
                sp.max.widths$se,sp.min.widths$se,sp.max.bowl.widths$se,sp.min.bowl.widths$se,sp.heights$se)
colnames(trait.se)<-c("x.cent","xs.disp","ys.cent","ys.disp",
                      "max.width","min.width","ss.max.width","ss.min.width","height")
trait.se[is.nan(trait.se)]<-NA
#write.csv(trait.dat,"trait.mat.csv")
#write.csv(trait.se,"trait.mat.se.csv")