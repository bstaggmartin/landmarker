##Data Wrangling Functions##

#Function for generating covariate matrix from raw data...
gen.covar.mat<-function(x,covar.start,lm.start,spec.lab="spec",leaf.lab="leaf"){
  EFN_att<-grep("EFN\\d+.type|EFN\\d+.side",colnames(x))
  splits<-strsplit(as.character(x[,spec.lab]),"_")
  x$sp<-transpose(splits)[[1]]
  x$ind<-transpose(splits)[[2]]
  x_covar<-x[seq(1,nrow(x)-2,3),c(ncol(x)-1,ncol(x),covar.start:(lm.start-1),EFN_att)]
  spec<-which(colnames(x_covar)==spec.lab)
  x_covar<-x_covar[,-spec]
  rownames(x_covar)<-x[,leaf.lab][seq(1,nrow(x)-2,3)]
  return(x_covar)
}
#Function for coercing raw data into helpful array format...
arrange.data<-function(x,lm.start,leaf.lab="leaf"){
  lms.cols<-lm.start:ncol(x)
  lms.cols<-lms.cols[!(lms.cols%in%grep("EFN\\d+.type|EFN\\d+.side",colnames(x)))]
  lms.rows<-rep(1:nrow(x),rep(c(1,1,0),nrow(x)/3))
  x.cleaned<-as.numeric(as.vector(as.matrix(t(x[lms.rows,lms.cols]))))
  x_sep<-array(data=x.cleaned,dim=c(length(lms.cols),2,nrow(x)/3),
               dimnames=list(colnames(x[,lms.cols]),c("x","y"),x[,leaf.lab][seq(1,nrow(x)-2,3)]))
  return(x_sep)
}
#Function for splitting each leaf into left and right half, reflecting the right half...
split.leaves<-function(x,x_covar,EFN.info.start){
  #left points and covariates...
  x.L<-x
  dimnames(x.L)[[3]]<-paste(dimnames(x)[[3]],".L",sep="")
  missing.Ls<-which(is.na(x["l.max","x",]))
  x.L<-x.L[,,-missing.Ls]
  L.EFN.inds<-which(x_covar=="L",arr.ind=T)
  L.EFN.nos<-sapply(strsplit(colnames(x_covar[L.EFN.inds[,2]]),"\\."),'[',1)
  L.EFN.nos<-substr(L.EFN.nos,4,nchar(L.EFN.nos))
  L.EFN.inds[,2]<-as.numeric(L.EFN.nos)
  L.EFN.inds<-L.EFN.inds[order(L.EFN.inds[,1]),]
  L.EFN.inds[,1]<-paste(rownames(x_covar)[L.EFN.inds[,1]],".L",sep="")
  x_covar.L<-x_covar
  rownames(x_covar.L)<-paste(rownames(x_covar),".L",sep="")
  x_covar.L<-x_covar.L[-missing.Ls,]
  for(i in dimnames(x.L)[[3]]){
    L.EFNs<-paste("EFN",L.EFN.inds[which(L.EFN.inds[,1]==i),2],"\\.",sep="")
    lms.to.keep<-c(1,2,3,grep(paste(L.EFNs,collapse="|"),dimnames(x)[[1]]))
    x.L[!(1:dim(x)[1]%in%lms.to.keep),,i]<-NA
    covars.to.keep<-c(1:(EFN.info.start-1),grep(paste(L.EFNs,collapse="|"),colnames(x_covar.L)))
    x_covar.L[i,!(1:ncol(x_covar.L)%in%covars.to.keep)]<-NA
  }
  x.L<-x.L[-4,,]
  dimnames(x.L)[[1]][3]<-"side.max"
  #right points and covariates...
  x.R<-x
  dimnames(x.R)[[3]]<-paste(dimnames(x)[[3]],".R",sep="")
  missing.Rs<-which(is.na(x["r.max","x",]))
  x.R<-x.R[,,-missing.Rs]
  R.EFN.inds<-which(x_covar=="R",arr.ind=T)
  R.EFN.nos<-sapply(strsplit(colnames(x_covar[R.EFN.inds[,2]]),"\\."),'[',1)
  R.EFN.nos<-substr(R.EFN.nos,4,nchar(R.EFN.nos))
  R.EFN.inds[,2]<-as.numeric(R.EFN.nos)
  R.EFN.inds<-R.EFN.inds[order(R.EFN.inds[,1]),]
  R.EFN.inds[,1]<-paste(rownames(x_covar)[R.EFN.inds[,1]],".R",sep="")
  x_covar.R<-x_covar
  rownames(x_covar.R)<-paste(rownames(x_covar),".R",sep="")
  x_covar.R<-x_covar.R[-missing.Rs,]
  for(i in dimnames(x.R)[[3]]){
    R.EFNs<-paste("EFN",R.EFN.inds[which(R.EFN.inds[,1]==i),2],"\\.",sep="")
    lms.to.keep<-c(1,2,4,grep(paste(R.EFNs,collapse="|"),dimnames(x)[[1]]))
    x.R[!(1:dim(x)[1]%in%lms.to.keep),,i]<-NA
    covars.to.keep<-c(1:(EFN.info.start-1),grep(paste(R.EFNs,collapse="|"),colnames(x_covar.R)))
    x_covar.R[i,!(1:ncol(x_covar.R)%in%covars.to.keep)]<-NA
  }
  x.R<-x.R[-3,,]
  dimnames(x.R)[[1]][3]<-"side.max"
  x.R[,"x",]<-x.R[,"x",]*-1
  #combine points and covariates with each half as separate specimens...
  new_x<-abind(x.L,x.R,along=3)
  new_x_covar<-rbind(x_covar.L,x_covar.R)
  return(list(x=new_x[,,order(as.numeric(sapply(strsplit(dimnames(new_x)[[3]],"\\."),'[',1)))],
              x_covar=new_x_covar[order(as.numeric(sapply(strsplit(rownames(new_x_covar),"\\."),'[',1))),]))
}
#Function for transforming each leaf to have the same location and rotation, with option to either standardized
#location by base lm (trans.by.cent=F) or centroid (trans.by.cent=T). Can optionally scale leaves to be the same
#size as well (scale=T/F)...
proc.fit<-function(x,trans.by.cent=F,scale=T){
  #define points that form basis for transformations...
  pts=c("tip","base","side.max")
  #translate...
  x_loc<-x
  x_loc[,"x",]<-x_loc[,"x",]-rep(apply(x[pts,"x",],2,mean),each=dim(x)[1])
  x_loc[,"y",]<-x_loc[,"y",]-rep(apply(x[pts,"y",],2,mean),each=dim(x)[1])
  #scale...
  scale_fac<-sqrt(apply(rbind(x_loc[pts,"x",],x_loc[pts,"y",])^2,2,mean)*2)
  x_scale<-x_loc/rep(scale_fac,each=2*dim(x)[1])
  #generate "ideal leaf" lm config...
  tips<-split(x_scale["tip",,],rep(1:dim(x_scale)[3],each=dim(x_scale)[2]))
  bases<-split(x_scale["base",,],rep(1:dim(x_scale)[3],each=dim(x_scale)[2]))
  side.maxes<-split(x_scale["side.max",,],rep(1:dim(x_scale)[3],each=dim(x_scale)[2]))
  length<-mean(mapply(dist.p,p1=tips,p2=bases))
  side.ys<-split(mapply(proj,p0=bases,p1=side.maxes,p2=tips),rep(1:dim(x_scale)[3],each=dim(x_scale)[2]))
  side.x<-mean(mapply(dist.p,p1=side.maxes,p2=side.ys))
  side.y<-mean(mapply(dist.p,p1=side.ys,p2=bases))
  ideal.config<-matrix(c(0,0,side.x,length,0,side.y),nrow=3)
  rownames(ideal.config)<-c("tip","base","side.max")
  colnames(ideal.config)<-c("x","y")
  ideal.config[,"x"]<-ideal.config[,"x"]-mean(ideal.config[,"x"])
  ideal.config[,"y"]<-ideal.config[,"y"]-mean(ideal.config[,"y"])
  #rotate...
  num<-apply(ideal.config[,"x"]*x_scale[1:3,"y",]-ideal.config[,"y"]*x_scale[1:3,"x",],2,sum)
  denom<-apply(ideal.config[,"y"]*x_scale[1:3,"y",]+ideal.config[,"x"]*x_scale[1:3,"x",],2,sum)
  vec<-ifelse(num<=0&denom<=0,"nn","pp")
  vec<-ifelse(num>=0&denom<=0,"np",vec)
  vec<-ifelse(num<=0&denom>=0,"pn",vec)
  rot_fac<-atan(num/denom)*-1
  rot_fac<-ifelse(vec=="nn",rot_fac+pi,rot_fac)
  rot_fac<-ifelse(vec=="np",rot_fac-pi,rot_fac)
  x_rot<-x_scale
  for(i in 1:dim(x_scale)[1]){
    landmarks<-split(x_scale[i,,],rep(1:dim(x_scale)[3],each=dim(x_scale)[2]))
    x_rot[i,,]<-mapply(rotate,p=landmarks,ang=rot_fac)
  }
  #optional final translate...
  if(trans.by.cent==F){
    x_rot[,"x",]<-x_rot[,"x",]-rep(x_rot["base","x",],each=dim(x)[1])
    x_rot[,"y",]<-x_rot[,"y",]-rep(x_rot["base","y",],each=dim(x)[1])
  }
  #optional final scale...
  if(scale==F){
    x_rot<-x_rot*rep(scale_fac,each=2*dim(x)[1])
  }
  return(x_rot)
}
#Function to extract lms (specified by name(s) or row number(s) in array--uses grep if provided w/ names!) and coerce
#them into list of coordinates with the name of each list element designated 'leafID landmarkID'...
extract.points<-function(x,p){
  if(!is.numeric(p)){
    p<-grep(p,dimnames(x)[[1]])
  }
  lst<-split(unlist(split(x[p,,],rep(1:length(p),dim(x)[2]))),rep(1:(dim(x)[3]*length(p)),each=2))
  names(lst)<-paste(rep(dimnames(x)[[3]],length(p)),rep(dimnames(x)[[1]][p],each=dim(x)[3]))
  return(lst)
}
#Function to extract x coordinates from list of coordinate pairs output by extract.points()...
extract.xs<-function(x){
  sapply(x,'[',1)
}
#Function to extract y coordinates form list of coordinate pairs output by extract.points()...
extract.ys<-function(x){
  sapply(x,'[',2)
}
#Function to create factor vector corresponding to a list of coordinates (one produced by either extract.xs/ys or 
#extract.points) for a given whole specimen variable--i.e., species or individuals...
spec.var<-function(x,x_covar,var){
  return(as.factor(x_covar[sapply(strsplit(names(x),' '),'[',1),var]))
}
#Function to create factor vector corresponding to a list of coordinates (one produced by either extract.xs/ys or 
#extract.points) for a given EFN variable--i.e., EFN types or sides (left vs. right)...
EFN.var<-function(x,x_covar,var){
  return(as.factor(x_covar[
    cbind(sapply(strsplit(names(x),' '),'[',1),
          paste(sapply(strsplit(sapply(strsplit(names(x),' '),'[',2),'\\.'),'[',1),'.',var,sep=''))]))
}
