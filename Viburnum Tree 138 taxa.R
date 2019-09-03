rm(list=ls())
def.par<-par(no.readonly=T)
library(phytools)

#cleaning up names to only say species by splitting tip labels based on under-score char (there was a
#missing underscore for V. stellato-tomentosum), then taking the second element of each split
Vib.Tree<-read.nexus("data/vib_tree_pollen_maxcred.txt")
Vib.Tree$tip.label[23]<-"Viburnum_stellato-tomentosum"
names.splits<-strsplit(Vib.Tree$tip.label,"_")
for(i in 1:length(names.splits)){
  Vib.Tree$tip.label[i]<-names.splits[[i]][2]
}
Vib.Tree<-force.ultrametric(Vib.Tree)
names(sp.heights$mean)[!(names(sp.heights$mean)%in%Vib.Tree$tip.label)]
Vib.Tree<-bind.tip(Vib.Tree,"hengshanicum",where=fastMRCA(Vib.Tree,"mullaha","foetidum"))
Vib.Tree<-bind.tip(Vib.Tree,"leiocarpum",where=fastMRCA(Vib.Tree,"inopinatum","sambucinum"))

pruned.tree<-drop.tip(Vib.Tree,
                      grep(paste(paste("^",c(names(sp.heights$mean)),"$",sep=""),collapse="|"),
                           Vib.Tree$tip.label,invert=T))
pruned.tree<-ladderize(pruned.tree)

pruned.tree$edge.length<-pruned.tree$edge.length/4

library(cairoDevice)
png(width=24*300,height=24*300,"test.png",bg="transparent",type="cairo")
seq<-seq(0.005,0.995,0.005)
par(xpd=T)
plot(pruned.tree,label.offset=18,cex=8,x.lim=50,edge.width=22.5,y.lim=50,align.tip.label=T)
test<-get("last_plot.phylo",envir=.PlotPhyloEnv)
offset=1
text(x=test$xx[1]+offset+base.pos*5,y=0.5,labels="base",cex=6,srt=90,adj=c(1,0.5))
lines(x=rep(test$xx[1]+offset+base.pos*5,2),y=c(0.75,max(test$yy)+1),lwd=4,lty=2)
text(x=test$xx[1]+offset+tip.pos*5,y=0.5,labels="tip",cex=6,srt=90,adj=c(1,0.5))
lines(x=rep(test$xx[1]+offset+tip.pos*5,2),y=c(0.75,max(test$yy)+1),lwd=4,lty=2)
for(i in 1:length(sp.EFNs.ys.shape1)){
  ys<-c(test$yy[order(pruned.tree$tip.label)][i],
        0.9*dbeta(seq,sp.EFNs.ys.shape1[i],sp.EFNs.ys.shape2[i])/
          max(dbeta(seq,sp.EFNs.ys.shape1[i],sp.EFNs.ys.shape2[i]))+
          test$yy[order(pruned.tree$tip.label)][i],
        test$yy[order(pruned.tree$tip.label)][i])-0.25
  xs<-c(offset+test$xx[1],5*seq+offset+test$xx[1],5*1+offset+test$xx[1])
  lines(ys~xs,lwd=2)
  polygon(ys[c((1:length(seq))+2,1:2)]~xs[c((1:length(seq))+2,1:2)],
          col=alter.cols(round(test$yy[order(pruned.tree$tip.label)][i]),0.25),lwd=3)
}
text(y=max(test$yy)+1.4,x=test$xx[1]+offset+0.3,cex=8,labels="Baso-apical Distribution",srt=30,adj=0)
offset=7
text(x=test$xx[1]+offset+midrib.pos*5,y=0.5,labels="midrib",cex=6,srt=90,adj=c(1,0.5))
lines(x=rep(test$xx[1]+offset+midrib.pos*5,2),y=c(0.75,max(test$yy)+1),lwd=4,lty=2)
text(x=test$xx[1]+offset+marg.pos*5,y=0.5,labels="margin",cex=6,srt=90,adj=c(1,0.5))
lines(x=rep(test$xx[1]+offset+marg.pos*5,2),y=c(0.75,max(test$yy)+1),lwd=4,lty=2)
for(i in 1:length(sp.EFNs.xs.shape1)){
  ys<-c(test$yy[order(pruned.tree$tip.label)][i],
        0.9*dbeta(seq,sp.EFNs.xs.shape1[i],sp.EFNs.xs.shape2[i])/
          max(dbeta(seq,sp.EFNs.xs.shape1[i],sp.EFNs.xs.shape2[i]))+
          test$yy[order(pruned.tree$tip.label)][i],
        test$yy[order(pruned.tree$tip.label)][i])-0.25
  xs<-c(offset+test$xx[1],5*seq+offset+test$xx[1],5*1+offset+test$xx[1])
  lines(ys~xs,lwd=2)
  polygon(ys[c((1:length(seq))+2,1:2)]~xs[c((1:length(seq))+2,1:2)],
          col=alter.cols(round(test$yy[order(pruned.tree$tip.label)][i]),0.25),lwd=3)
}
text(y=max(test$yy)+1.4,x=test$xx[1]+offset+0.3,cex=8,labels="Medio-lateral Distribution",srt=30,adj=0)
offset=13
for(i in 1:length(sp.max.widths$mean)){
  draw.ellipse(x=test$xx[1]+0.5+offset,y=test$yy[order(pruned.tree$tip.label)][i],
               a=exp(sp.max.widths$mean[i])/2,b=exp(sp.min.widths$mean[i])/2,angle=0,
               lwd=3,col="chartreuse3")
  draw.ellipse(x=test$xx[1]+0.5+offset,y=test$yy[order(pruned.tree$tip.label)][i],
               a=exp(sp.max.bowl.widths$mean[i])/2,b=exp(sp.min.bowl.widths$mean[i])/2,angle=0,
               lwd=3,col="darkorange3")
}
offset=14.5
lines(y=c(0,0),x=test$xx[1]+offset+c(-0.5,0.5),lwd=15)
text(x=test$xx[1]+offset,y=-0.5,adj=0.5,labels="1 mm",cex=6)
text(y=max(test$yy)+1.4,x=test$xx[1]+offset-0.3,cex=8,labels="Shape and Size",srt=30,adj=0)
offset=15
for(i in 1:length(sp.max.widths$mean)){
  polygon(x=c(test$xx[1]-exp(sp.max.widths$mean)[i]/2,test$xx[1]-exp(sp.max.bowl.widths$mean)[i]/2,
              test$xx[1]+exp(sp.max.bowl.widths$mean)[i]/2,test$xx[1]+exp(sp.max.widths$mean)[i]/2)+offset+0.5,
          y=c(test$yy[order(pruned.tree$tip.label)][i]-exp(sp.heights$mean[i])/2,
              test$yy[order(pruned.tree$tip.label)][i]+exp(sp.heights$mean[i])/2,
              test$yy[order(pruned.tree$tip.label)][i]+exp(sp.heights$mean[i])/2,
              test$yy[order(pruned.tree$tip.label)][i]-exp(sp.heights$mean[i])/2),
          lwd=3,col="chartreuse3")
  draw.ellipse(x=test$xx[1]+offset+0.5,y=test$yy[order(pruned.tree$tip.label)][i]+exp(sp.heights$mean[i])/2,
               a=exp(sp.max.bowl.widths$mean[i])/2,b=exp(sp.max.bowl.widths$mean[i])/10,angle=0,
               lwd=3,col="darkorange3")
}
offset=17
text(y=max(test$yy)+1.4,x=test$xx[1]+offset,cex=8,labels="Latitude",srt=30,adj=0)
#points(x=rep(test$xx[1]+offset,length(sp.seasonality)),y=test$yy[order(pruned.tree$tip.label)],
       #pch=16,col=rgb(0,0.75-0.5*std.sp.seasonality,0+0.75*std.sp.seasonality),cex=12)
#points(y=c(0.2,0.2),x=c(-0.5,0.5)+test$xx[1]+offset,col=rgb(0,c(0.75,0),c(0.25,0.75)),pch=16,cex=12)
points(x=rep(test$xx[1]+offset,length(sp.temp2)),y=test$yy[order(pruned.tree$tip.label)],
       pch=16,col=rgb(rbind(c(0,0.25,0.75),c(0,0.75,0)))[sp.temp2+1],cex=12)
points(y=c(0.2,0.2),x=c(-0.5,0.5)+test$xx[1]+offset,col=rgb(0,c(0.25,0.75),c(0.75,0)),pch=16,cex=12)
#text(x=test$xx[1]+offset+c(-0.5,0.5),y=c(-0.3,-0.3),adj=1,labels=c("Aseasonal","Seasonal"),cex=4,srt=30)
text(x=test$xx[1]+offset+c(-0.75,0.75),y=c(-0.5,-0.5),adj=0.5,labels=c(">23.5","<23.5"),cex=6)
dev.off()



par(mfrow=c(2,2))
plotTree.wBars(pruned.tree,sp.EFNs.xs.shape1)
plotTree.wBars(pruned.tree,sp.EFNs.xs.shape2)
plotTree.wBars(pruned.tree,sp.EFNs.ys.shape1)
plotTree.wBars(pruned.tree,sp.EFNs.ys.shape2)
par(def.par)

phenogram(pruned.tree,exp(sp.means)[names(sp.means)%in%pruned.tree$tip.label])
phenogram(pruned.tree,sp.EFNs.nos[names(sp.EFNs.nos)%in%pruned.tree$tip.label])

seq<-seq(0.01,0.99,0.01)
y.adv<-2
x.adv<-0.1
plot(c(0,dbeta(seq,sp.EFNs.xs.shape1[1],sp.EFNs.xs.shape2[1]),0)~c(0,seq,1),
     xlim=c(0,1.5+length(sp.EFNs.xs.shape1)*x.adv),ylim=c(0,10+length(sp.EFNs.xs.shape1)*y.adv),type="l",
     xaxt='n',xlab="Medio-lateral EFN Position",yaxt='n',ylab="Probability Density",line=1)
axis(1,c(midrib.pos,marg.pos),c("midrib","margin"))
polygon(c(dbeta(seq,sp.EFNs.xs.shape1[1],sp.EFNs.xs.shape2[1]),0,0)~c(seq,1,0),col=alter.cols(1,0.25),border=NA)
text(1,0,names(sp.EFNs.xs.shape1[1]),adj=-0.1,cex=0.6)
for(i in 2:length(sp.EFNs.xs.shape1)){
  ys<-dbeta(seq,sp.EFNs.xs.shape1[i],sp.EFNs.xs.shape2[i])+((i-1)*y.adv)
  xs<-seq+((i-1)*x.adv)
  lines(c((i-1)*y.adv,ys,(i-1)*y.adv)~c((i-1)*x.adv,xs,(i-1)*x.adv+1))
  polygon(c(ys,(i-1)*y.adv,(i-1)*y.adv)~c(xs,(i-1)*x.adv+1,(i-1)*x.adv),col=alter.cols(i,0.25),border=NA)
  text((i-1)*x.adv+1,(i-1)*y.adv,names(sp.EFNs.xs.shape1[i]),adj=-0.1,cex=0.6)
}


plot(c(0,dbeta(seq,sp.EFNs.ys.shape1[1],sp.EFNs.ys.shape2[1]),0)~c(0,seq,1),
     xlim=c(0,1.5+length(sp.EFNs.ys.shape1)*x.adv),ylim=c(0,10+length(sp.EFNs.ys.shape1)*y.adv),type="l",
     xaxt='n',xlab="Baso-apical EFN Position",yaxt='n',ylab="Probability Density",line=1)
axis(1,c(base.pos,tip.pos),c("base","tip"))
polygon(c(dbeta(seq,sp.EFNs.ys.shape1[1],sp.EFNs.ys.shape2[1]),0,0)~c(seq,1,0),col=alter.cols(1,0.25),border=NA)
text(1,0,names(sp.EFNs.xs.shape1[1]),adj=-0.1,cex=0.6)
for(i in 2:length(sp.EFNs.ys.shape1)){
  ys<-dbeta(seq,sp.EFNs.ys.shape1[i],sp.EFNs.ys.shape2[i])+((i-1)*y.adv)
  xs<-seq+((i-1)*x.adv)
  lines(c((i-1)*y.adv,ys,(i-1)*y.adv)~c((i-1)*x.adv,xs,(i-1)*x.adv+1))
  polygon(c(ys,(i-1)*y.adv,(i-1)*y.adv)~c(xs,(i-1)*x.adv+1,(i-1)*x.adv),col=alter.cols(i,0.25),border=NA)
  text((i-1)*x.adv+1,(i-1)*y.adv,names(sp.EFNs.xs.shape1[i]),adj=-0.1,cex=0.6)
}


plot(pruned.tree)
std.sp.lats<-sp.lats-min(sp.lats);std.sp.lats<-std.sp.lats/max(std.sp.lats)
points(x=rep(test$xx[1]+10,length(sp.lats)),y=test$yy[order(pruned.tree$tip.label)],
       pch=16,col=rgb(0.1,1-0.5*std.sp.lats,0.5+0.5*std.sp.lats))

plot(pruned.tree)
std.sp.seasonality<-sp.seasonality-min(sp.seasonality);std.sp.seasonality<-std.sp.seasonality/max(std.sp.seasonality)
points(x=rep(test$xx[1]+10,length(sp.seasonality)),y=test$yy[order(pruned.tree$tip.label)],
       pch=16,col=rgb(0.1,1-0.5*std.sp.seasonality,0.5+0.5*std.sp.seasonality))

sp.temp<-ifelse(sp.lats<25,0,1)

plot(pruned.tree)
test<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(x=rep(test$xx[1]+5,length(sp.seasonality)),y=test$yy[order(pruned.tree$tip.label)],
       pch=16,col=rgb(0.1,0.75-0.5*sp.temp,0+0.75*sp.temp))

plot(0,xlim=c(0,15),ylim=c(0,10),col="white",asp=1)
for(i in 1:n.spp){
  xx<-std.sp.lats[i]*15
  yy<-std.sp.seasonality[i]*10
  polygon(x=c(xx-exp(sp.max.widths$mean)[i]/2,xx-exp(sp.max.bowl.widths$mean)[i]/2,
              xx+exp(sp.max.bowl.widths$mean)[i]/2,xx+exp(sp.max.widths$mean)[i]/2),
          y=c(yy,yy+exp(sp.heights$mean)[i],yy+exp(sp.heights$mean)[i],yy),
          col=alter.cols("chartreuse3",0.6,name=T),lwd=1)
  draw.ellipse(x=xx,y=yy+exp(sp.heights$mean)[i],
               a=exp(sp.max.bowl.widths$mean)[i]/2,b=0.2*exp(sp.heights$mean)[i],angle=0,
               col=alter.cols("darkorange3",0.6,name=T),lwd=1)
}

plot(0,xlim=c(0,15),ylim=c(0,10),col="white",asp=1)
draw.ellipse(x=std.sp.lats*15,y=std.sp.seasonality*10,
             a=exp(sp.max.widths$mean)/2,b=exp(sp.min.widths$mean)/2,angle=rep(0,n.spp),
             col=alter.cols("chartreuse3",0.6,name=T),lwd=1)
draw.ellipse(x=std.sp.lats*15,y=std.sp.seasonality*10,
             a=exp(sp.max.bowl.widths$mean)/2,b=exp(sp.min.bowl.widths$mean)/2,angle=rep(0,n.spp),
             col=alter.cols("darkorange3",0.6,name=T),lwd=1)
