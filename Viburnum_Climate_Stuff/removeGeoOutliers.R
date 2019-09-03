removeGeoOutliers <- function(coords){
  nas<-which(is.na(coords$lon)|is.na(coords$lat))
  coords.no.nas<-coords[-nas,]
  if(nrow(coords.no.nas)>0){
    coords.no.nas$lon <- round(coords.no.nas$lon,4)
    coords.no.nas$lat <- round(coords.no.nas$lat,4)
    duplicates <- which(duplicated(coords.no.nas[, c('lon', 'lat')]))
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    wholenums <- as.numeric(rownames(subset(coords.no.nas, (is.wholenumber(lon*10) | is.wholenumber(lat*10)))))
    lon.out <- coords.no.nas[coords.no.nas$lon %in% unlist(boxplot.stats(coords.no.nas$lon, coef = 1.5, do.conf = TRUE, do.out = TRUE)$out), ]
    lat.out <- coords.no.nas[coords.no.nas$lat %in% unlist(boxplot.stats(coords.no.nas$lat, coef = 1.5, do.conf = TRUE, do.out = TRUE)$out), ]
    problem.coords<-c(duplicates,wholenums,lon.out,lat.out,nas)
    problem.coords<-problem.coords[!duplicated(problem.coords)]
  }else{
    problem.coords<-nas
  }
  for(i in 1:nrow(coords)){
    if(i%in%problem.coords){
      coords[i,"lat"]=NA
      coords[i,"lon"]=NA
    }else{
      coords[i,"lat"]<-round(coords[i,"lat"],4)
      coords[i,"lon"]<-round(coords[i,"lon"],4)
    }
  }
  return(coords)
}
