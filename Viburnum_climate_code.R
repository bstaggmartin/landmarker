rm(list=ls())

####################################
# 1. load packages ####
####################################
require(ape)
require(rgbif)
require(maps)
require(dismo)
require(ape)
require(rgeos)
require(maptools)
require(rgdal)
source('Viburnum_Climate_Stuff/removeGeoOutliers.R')

species<-Vib.Tree$tip.label
genera<-rep("Viburnum",length(species))
####################################
# 2 Function to download and clean coordinates ####
####################################


getCoords <- function(genus, epithet){
  # genus='Sidal'; epithet='hirtipes'
  # genus<-genera[i]; epithet<-species[i]
  cat('\n\n\n')
  
  xy <- gbif(genus, paste(epithet, '*', sep = ''), ntries = 10)
  #xy2 <- getConsortium(paste(genus, epithet))
  if(!is.null(xy)) xy <- xy[ , c('lat', 'lon')]
  
  xy <- removeGeoOutliers(xy)
  if(nrow(xy)>0){
    return(xy)
  }else{
    return(data.frame(lat=NA,lon=NA))
  }
}



####################################
# 3 Run the function ####
####################################


coords <- list()


for (i in (1:length(genera))[-which(species=="hengshanicum")]){
  species.temp<-paste(genera[i], species[i])
  cat(i,species.temp)
  if(!species.temp %in% names(list)){
    coords[[i]]<- data.frame(sp = species.temp, 
                             getCoords(genera[i], species[i]))
    names(coords)[i]<-paste(genera[i], species[i])
  }
}




#remove NAs


conds<-unlist(lapply(coords, function(x){all(is.na(x[,2])|is.na(x[,3]))}))
missingTaxa<-names(conds)[conds] #identify problematic taxa
missingTaxa

coords1<- coords[!conds] #remove problematic taxa
coordsGBIF <- do.call(rbind, coords1)
head(coordsGBIF)
row.names(coordsGBIF)<-NULL

####################################
# 4 make plots of GBIF data ####
####################################

coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum acerifolium"&
                           (coordsGBIF[,"lon"]> -45|coordsGBIF[,"lon"]< -110)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum awabuki"&coordsGBIF[,"lon"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum betulifolium"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum blandum"&coordsGBIF[,"lat"]>24),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum bracteatum"&coordsGBIF[,"lon"]> -15),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum carlesii"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum carlesii"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum chingii"&coordsGBIF[,"lon"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum caudatum"&coordsGBIF[,"lon"]> -50),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum cinnamomifolium"&coordsGBIF[,"lon"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum coriaceum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum cylindricum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum davidii"&
                           (coordsGBIF[,"lon"]>150|coordsGBIF[,"lon"]<75)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum dentatum"&
                           (coordsGBIF[,"lon"]< -105|coordsGBIF[,"lon"]> -50|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum dilatatum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum erosum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum farreri"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum foetidum"&coordsGBIF[,"lat"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum furcatum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum hallii"&coordsGBIF[,"lat"]< -24),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum hupehense"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum integrifolium"&coordsGBIF[,"lon"]<75),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum japonicum"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum kansuense"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum lantana"&coordsGBIF[,"lon"]< -15),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum lantanoides"&
                           (coordsGBIF[,"lon"]> -50|coordsGBIF[,"lon"]< -105)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum lautum"&coordsGBIF[,"lon"]> -50),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum lentago"&
                           (coordsGBIF[,"lon"]> -50|coordsGBIF[,"lon"]< -115)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum lobophyllum"&coordsGBIF[,"lat"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum luzonicum"&
                           (coordsGBIF[,"lon"]< -50|coordsGBIF[,"lat"]< -15)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum macrocephalum"&
                           (coordsGBIF[,"lat"]<0|coordsGBIF[,"lon"]<60)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum melanocarpum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum molle"&coordsGBIF[,"lon"]> -80),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum mongolicum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum mullaha"&coordsGBIF[,"lon"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum nudum"&
                           (coordsGBIF[,"lon"]< -105|coordsGBIF[,"lon"]> -50|coordsGBIF[,"lat"]<15)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum odoratissimum"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum opulus"&
                           (coordsGBIF[,"lon"]>60|coordsGBIF[,"lat"]<0|coordsGBIF[,"lon"]< -15)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum plicatum"&
                           (coordsGBIF[,"lon"]<90|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum propinquum"&coordsGBIF[,"lat"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum prunifolium"&
                           (coordsGBIF[,"lon"]> -50|coordsGBIF[,"lon"]< -115)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum rafinesquianum"&coordsGBIF[,"lat"]<35),]
#^restricting to N variety for now (sample from MN anyways...)
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum rhytidophyllum"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum rigidum"&
                           (coordsGBIF[,"lat"]>35|coordsGBIF[,"lat"]<20)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum rufidulum"&
                           (coordsGBIF[,"lon"]> -75|coordsGBIF[,"lon"]< -110)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum sargentii"&coordsGBIF[,"lon"]<90),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum schensianum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum setigerum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum sieboldii"&coordsGBIF[,"lon"]<90),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum subalpinum"&coordsGBIF[,"lon"]<60),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum suspensum"&
                           (coordsGBIF[,"lon"]<60|coordsGBIF[,"lat"]<0)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum tinus"&
                           (coordsGBIF[,"lon"]>60|coordsGBIF[,"lon"]< -30)),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum trilobum"&coordsGBIF[,"lat"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum veitchii"&coordsGBIF[,"lat"]<0),]
coordsGBIF<-coordsGBIF[!(coordsGBIF[,"sp"]=="Viburnum wrightii"&coordsGBIF[,"lon"]<60),]

coordsGBIF[,"sp"]<-as.character(coordsGBIF[,"sp"])
coordsGBIF<-rbind(coordsGBIF,c("Viburnum hengshanicum",27.3,112.8),c("Viburnum hengshanicum",30.2,118.2))
coordsGBIF[,"sp"]<-as.factor(coordsGBIF[,"sp"])
sp<-sort(unique(coordsGBIF$sp))

for (i in 1:length(sp)){
  xy<-coordsGBIF[coordsGBIF$sp==sp[i],]
  filenameTemp<-paste("Viburnum_Climate_Stuff/1_MapPlots/",
                      unlist(strsplit(as.character(sp[i]), split= ' '))[2],".png",sep="")
  png(filename=filenameTemp)  
  map("world", fill=TRUE,col="grey")#,ylim=c(-60,30),xlim=c(-110,-35))
  points(xy$lon,xy$lat,pch=16,col=2)
  text(-125,-50,as.character(sp[i]))
  text(-125,-60,paste(as.character(nrow(xy)),"collections"))
  dev.off()
}


####################################
# 5 Export GBIF Data ####
####################################

write.csv(coordsGBIF, file = 'Viburnum_Climate_Stuff/Viburnum_coords_gbif.csv', row.names = F)


####################################
# 5 Make ranges ####
####################################

# get climate means for species ## 
viburnum.coords<-read.csv('Viburnum_Climate_Stuff/Viburnum_coords_gbif.csv')
head(viburnum.coords)

# 01. Read in data

climate <- getData('worldclim', var='bio', res=2.5)


# make raster plots...
plot(climate[[4]],
     main="Seasonality",xlab="Longitude",ylab="Latitude",
     cex.axis=1.3,cex.lab=1.4,cex.main=1.5,col=terrain.colors(50))


# 02. Extract climate values
coords<-data.frame(lon=viburnum.coords$lon, lat=viburnum.coords$lat)
coordinates(coords)<-c("lon","lat")
Vib.climVals <- extract(climate, coords)
Vib.climVals <- data.frame(viburnum.coords, Vib.climVals)
head(Vib.climVals)

write.csv(Vib.climVals,'data/Viburnum.climateCoords.csv')

# 03. Get climate means by species
dat <- aggregate(Vib.climVals, list(Vib.climVals$sp), mean, na.rm = T)
# 04. Expor
write.csv(dat,'data/Viburnum.climateMeans.csv')
