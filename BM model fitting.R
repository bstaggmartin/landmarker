library(mvMORPH)

tmat<-as.matrix(read.csv("trait.mat.csv",row.names=1))
tsemat<-as.matrix(read.csv("trait.mat.se.csv",row.names=1))

std.tmat<-apply(tmat,2,function(x) {(x-mean(x,na.rm=T))/sd(x,na.rm=T)})
std.tsemat<-tsemat
for(i in 1:ncol(tsemat)){
  std.tsemat[,i]<-tsemat[,i]/sd(tmat[,i])
}


fit1<-mvBM(pruned.tree,std.tmat[,-9],std.tsemat[,-9],model="BM1",method="rpf",scale.height=T)
fit2<-mvBM(pruned.tree,std.tmat[,-9],std.tsemat[,-9],model="BM1",method="rpf",scale.height=T,
           params<-list(constraint="diagonal"))
fit3<-mvBM(pruned.tree,std.tmat[,-9],std.tsemat[,-9],model="BM1",method="rpf",scale.height=T,
           params<-list(constraint="equal"))
LRT(fit1,fit3)

heatmap(fit1$sigma)

red.fit<-mvBM(pruned.tree,std.tmat[,1:4],std.tsemat[,1:4],model="BM1",method="rpf",scale.height=T)
red.fit2<-mvBM(pruned.tree,std.tmat[,1:4],std.tsemat[,1:4],model="BM1",method="rpf",scale.height=T,
               params<-list(constraint="diagonal"))
best.l.fit<-mvBM(pruned.tree,std.tmat[,1:4],std.tsemat[,1:4],model="BM1",method="rpf",scale.height=T,
               params<-list(constraint="equal"))
LRT(red.fit,red.fit3)

sred.fit1<-mvBM(pruned.tree,std.tmat[,5:8],std.tsemat[,5:8],model="BM1",method="rpf",scale.height=T)
sred.fit2<-mvBM(pruned.tree,std.tmat[,5:8],std.tsemat[,5:8],model="BM1",method="rpf",scale.height=T,
               params<-list(constraint="diagonal"))
sred.fit3<-mvBM(pruned.tree,std.tmat[,5:8],std.tsemat[,5:8],model="BM1",method="rpf",scale.height=T,
               params<-list(constraint="equal"))
LRT(sred.fit1,sred.fit3)

loc.map<-make.simmap(tree=pruned.tree,x=sp.temp,nsim=100,model='ARD')
plotSimmap(loc.map)

cols<-setNames(c("green","blue"),c(0,1))
obj<-summary(loc.map)
plot(obj,colors=cols,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=10,y=10,vertical=FALSE,shape="circle")

best.s.fit<-mvOU(pruned.tree,std.tmat[,5:8],std.tsemat[,5:8],model="OU1",method="rpf",scale.height=T)
          #params<-list(constraint="equal"))
#support for non-equal model
mvOU(pruned.tree,std.tmat[,1:4],std.tsemat[,1:4],model="OU1",method="rpf",scale.height=T)
     #params<-list(constraint="equal"))
#bad support whether equal or not
mvEB(pruned.tree,std.tmat[,5:8],std.tsemat[,5:8],method="rpf",scale.height=T)
#pretty poor support...

#so, to summarize: EFN locations and size/shape seem to evolve independently (neat) yet are highly  correlated amongst
#themselves, and an OU best seems to explain evolution of EFN size and a BM evolution of EFN locations
theta<-exp((best.s.fit$theta)*apply(tmat[,5:8],2,sd,na.rm=T)+apply(tmat[,5:8],2,mean,na.rm=T))

#Q<-matrix(c(-0.1,0.1,-0.01,0.01),byrow=T,ncol=2);rownames(Q)<-c(0,1);colnames(Q)<-c(0,1)
map<-make.simmap(pruned.tree,x=sp.temp2,model='ARD',nsim=1)
mvBM(map,std.tmat[,5:8],std.tsemat[,5:8],model="BMM",method="rpf",scale.height=T,
     params<-list(constraint="equal"))
#this is best balance between simplicity and fit
mvBM(map,std.tmat[,5:8],std.tsemat[,5:8],model="BMM",method="rpf",scale.height=T)


mvBM(map,std.tmat[,1:4],std.tsemat[,1:4],model="BMM",method="rpf",scale.height=T,
     params<-list(constraint="equal"))
mvBM(map,std.tmat[,1:4],std.tsemat[,1:4],model="BMM",method="rpf",scale.height=T)
#this fits best...
mvOU(map,std.tmat[,1:4],std.tsemat[,1:4],model="OUM",method="rpf",scale.height=T)

##very preliminary model fitting, but results so far suggest that EFN location and size show opposing patterns
##EFN size evolution slows down in tropics, while locational evolution increases--particularly in terms of baso-apical
##distribution

##okay, so size conclusion disappears if you raise the tropicality  cutoff by just 1.5 degs...
##BUT location conclusion still passes with flying colors!

Vibminsplit<-minSplit(pruned.tree,testy$mcmc[50:501,c("node","bp")],method="sum")

#a quick rjMCMC seems to indicate no noticeable rate shifts in size evolution, but a dsitinct one for the base of
#lobata + succotinus, indicating a significant slow down in location evolution

sig<-fit1$sigma
colnames(sig)<-c("ML centroid","ML dispersion","BA centroid","BA dispersion",
                 "max width","min width","max ss width","min ss width")
rownames(sig)<-c("ML centroid","ML dispersion","BA centroid","BA dispersion",
                 "max width","min width","max ss width","min ss width")
sig<-sig[rev(c(6,5,8,7,4,3,2,1)),rev(c(1,2,3,4,7,8,5,6))]
png(width=2000,height=2000,file="heatmap.png",type='cairo',bg='transparent')
levelplot(sig,xaxt='n',yaxt='n',xlab='',ylab='',at=seq(-2,6),
          col.regions=rgb(colorRamp(c("deepskyblue3","green2","yellow","red"))(seq(0,1,length.out=100))/255),
          main="",las=2,scales=list(draw=FALSE),colorkey=list(width=15,axis.text=list(cex=8),space="left"))
dev.off()


fit1.2<-mvBM(pruned.tree,tmat[,-9],tsemat[,-9],model="BM1",method="rpf",scale.height=T)
fit3.2<-mvBM(pruned.tree,tmat[,-9],tsemat[,-9],model="BM1",method="rpf",scale.height=T,
           params<-list(constraint="equal"))
LRT(fit1.2,fit3.2)

fit1.3<-mvBM(map,tmat[,-9],tsemat[,-9],model="BMM",method="rpf",scale.height=T)
fit1.4<-mvBM(map,std.tmat[,-9],std.tsemat[,-9],model="BMM",method="rpf",scale.height=T)


fit1.5<-mvBM(map,std.tmat[,-9],std.tsemat[,-9],model="BMM",method="rpf",scale.height=T,
             param=list(constraint="shared"))

user_const<-array(rep(17),dim=c(8,8,2))
user_const[cbind(1:8,1:8,rep(c(1,2),each=8))]<-1:16

fit1.6<-mvBM(map,std.tmat[,-9],std.tsemat[,-9],model="BMM",method="rpf",scale.height=T,
             param=list(constraint=user_const))

set.seed(123)
map<-make.simmap(pruned.tree,x=sp.temp2,model='ARD',nsim=100)
dist.BMM.fit<-list()
dist.BM1.fit<-list()
s.BMM.fit<-list()
s.BM1.fit<-list()
for(i in 1:100){
  dist.BMM.fit[[i]]<-mvBM(map[[i]],std.tmat[,1:4],std.tsemat[,1:4],model="BMM",method="rpf",scale.height=T,echo=F)
  dist.BM1.fit[[i]]<-mvBM(pruned.tree,std.tmat[,1:4],std.tsemat[,1:4],model="BM1",method="rpf",scale.height=T,echo=F)
  s.BMM.fit[[i]]<-mvBM(map[[i]],std.tmat[,5:8],std.tsemat[,5:8],model="BMM",method="rpf",scale.height=T,echo=F)
  s.BM1.fit[[i]]<-mvBM(pruned.tree,std.tmat[,5:8],std.tsemat[,5:8],model="BM1",method="rpf",scale.height=T,echo=F)
  print(i)
}
#saveRDS(dist.BMM.fit,"dist.BMM.fit.obj")
#saveRDS(dist.BM1.fit,"dist.BM1.fit.obj")
#saveRDS(s.BMM.fit,"s.BMM.fit.obj")
#saveRDS(s.BM1.fit,"s.BM1.fit.obj")
dist.BMM.fit<-readRDS("dist.BMM.fit.obj")
dist.BM1.fit<-readRDS("dist.BM1.fit.obj")
s.BMM.fit<-readRDS("s.BMM.fit.obj")
s.BM1.fit<-readRDS("s.BM1.fit.obj")

mean(unlist(sapply(dist.BMM.fit,'[','LogLik'))-unlist(sapply(dist.BM1.fit,'[','LogLik')))
hist(unlist(sapply(s.BM1.fit,'[','AICc'))-unlist(sapply(s.BMM.fit,'[','AICc')))

plot(sort(unlist(mapply(LRT,model1=dist.BMM.fit,model2=dist.BM1.fit)["pval",])))
abline(h=0.05,col="yellow")
abline(h=0.01,col='orange')
abline(h=0.001,col='red')

plot(sort(unlist(mapply(LRT,model1=s.BMM.fit,model2=s.BM1.fit)["pval",])))
abline(h=0.05,col="yellow")
abline(h=0.01,col='orange')
abline(h=0.001,col='red')

s.sigs<-sapply(s.BMM.fit,'[','sigma')
dist.sigs<-sapply(dist.BMM.fit,'[','sigma')
s.sig.difs<-list()
for(i in 1:length(sigs)){
  s.sig.difs[[i]]<-s.sigs[[i]][,,2]-sigs[[i]][,,1]
}
rate.dif<-unlist(sapply(sig.difs,'[',c(4,4)))
hist(rate.dif)

s.temp.sigs<-array(NA,dim=c(4,4,100))
s.trop.sigs<-array(NA,dim=c(4,4,100))
dist.temp.sigs<-array(NA,dim=c(4,4,100))
dist.trop.sigs<-array(NA,dim=c(4,4,100))
for(i in 1:length(s.sigs)){
  s.temp.sigs[,,i]<-s.sigs[[i]][,,1]
  s.trop.sigs[,,i]<-s.sigs[[i]][,,2]
  dist.temp.sigs[,,i]<-dist.sigs[[i]][,,1]
  dist.trop.sigs[,,i]<-dist.sigs[[i]][,,2]
}

avg.temp.dist.rates<-apply(dist.temp.sigs,c(1,2),mean)
avg.temp.dist.rates.se<-apply(dist.temp.sigs,c(1,2),sd)
avg.temp.s.rates<-apply(s.temp.sigs,c(1,2),mean)
avg.temp.s.rates.se<-apply(s.temp.sigs,c(1,2),sd)
avg.trop.dist.rates<-apply(dist.trop.sigs,c(1,2),mean)
avg.trop.dist.rates.se<-apply(dist.trop.sigs,c(1,2),sd)
avg.trop.s.rates<-apply(s.trop.sigs,c(1,2),mean)
avg.trop.s.rates.se<-apply(s.trop.sigs,c(1,2),sd)

s.rate.difs<-s.trop.sigs-s.temp.sigs
png(width=2300,height=2000,file="s.rates.png",bg='transparent',type='cairo')
par(mar=c(30,30,10,10),xpd=T)
plot(0,col="white",xlim=c(-14,8),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
axis(1,cex.axis=8,line=5,col=rgb(0,0,0,0),col.ticks=rgb(0,0,0,0))
axis(2,cex.axis=8,line=3,col=rgb(0,0,0,0),col.ticks=rgb(0,0,0,0))
mtext(side=1,line=16,cex=11,text=bquote(sigma[tropical] - sigma[temperate]))
mtext(side=2,line=13,cex=8,text="Bootstrapped Probability Density")
box(lwd=12)
lines(c(0,1)~c(0,0),lwd=12,lty=2)
polygon(density(s.rate.difs[1,1,],adjust=3),col=alter.cols(1,0.25),lwd=6)
polygon(density(s.rate.difs[2,2,],adjust=3),col=alter.cols(3,0.25),lwd=6)
polygon(density(s.rate.difs[3,3,],adjust=3),col=alter.cols(5,0.25),lwd=6)
polygon(density(s.rate.difs[4,4,],adjust=3),col=alter.cols(7,0.25),lwd=6)
legend("topleft",pt.bg=alter.cols(c(1,3,5,7),0.25),pch=22,bty='n',pt.cex=8,cex=6,pt.lwd=6,
       legend=c("max width","min width","max secretory surface width","min secretory surface width"))
dev.off()

dist.rate.difs<-dist.trop.sigs-dist.temp.sigs
png(width=2300,height=2000,file="dist.rates.png",bg='transparent',type='cairo')
par(mar=c(30,30,10,10),xpd=T)
plot(0,col="white",xlim=c(-2,20),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='')
axis(1,cex.axis=8,line=5,col=rgb(0,0,0,0),col.ticks=rgb(0,0,0,0))
axis(2,cex.axis=8,line=3,col=rgb(0,0,0,0),col.ticks=rgb(0,0,0,0))
mtext(side=1,line=16,cex=11,text=bquote(sigma[tropical] - sigma[temperate]))
mtext(side=2,line=13,cex=8,text="Bootstrapped Probability Density")
box(lwd=12)
lines(c(0,1)~c(0,0),lwd=12,lty=2)
polygon(density(dist.rate.difs[1,1,],adjust=3),col=alter.cols(4,0.25),lwd=6)
polygon(density(dist.rate.difs[2,2,],adjust=3),col=alter.cols(6,0.25),lwd=6)
polygon(density(dist.rate.difs[3,3,],adjust=3),col=alter.cols(8,0.25),lwd=6)
polygon(density(dist.rate.difs[4,4,],adjust=3),col=alter.cols(10,0.25),lwd=6)
legend("topright",pt.bg=alter.cols(c(4,6,8,10),0.25),pch=22,bty='n',pt.cex=8,cex=6,pt.lwd=6,
       legend=c("margino-lateral centroid","margino-lateral dispersion",
                "baso-apical centroid","baso-apical dispersion"))
dev.off()
