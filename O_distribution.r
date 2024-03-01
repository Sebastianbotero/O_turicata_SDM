library(sf)
library(terra)
library(spatstat)
library(raster)
library(adehabitatHR)
library(XML)
library(rgdal)
library(gdalUtils)
library(ENMTools)

## Equal area projection
Albers<-st_crs("ESRI:102003")


### Background sampling points - Tick records through study area
bkg<-read.csv("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Ixodida.csv")
bkg1<- st_as_sf(bkg[!is.na(bkg$decimalLatitude),], coords = c("decimalLongitude","decimalLatitude"))
st_crs(bkg1)<-st_crs("EPSG:4326")
bkg<-st_transform(bkg1,Albers)


#### Area Of interest 1
AOI_wgs84<-read_sf("C:/Users/seboc/Box/Argasid_tick_project/SDM_GIS_Feb2024/AOI_EasternPop.shp")
AOI<-st_transform(AOI_wgs84,Albers)


#### Import Bios+ from CHELSA - include 19 bios plus near-surface relative humidity (hurs) and site water balance (swb)

setwd("C:/Users/seboc/Box/Argasid_tick_project/SDM_GIS_Feb2024/CHELSA_BIOS/")
bios<-list.files("C:/Users/seboc/Box/Argasid_tick_project/SDM_GIS_Feb2024/CHELSA_BIOS/")

bios <- lapply(bios[-c(10,11,18,19)], rast)
bios<-rast(bios)

bios1<-crop(bios,AOI_wgs84)
bios1<-terra::project(bios1,"ESRI:102003")
bios2<-mask(bios1,AOI)
plot(bios2[[6]])

a<-names(bios2)
a1<-strsplit(a,"_")
a2<-lapply(a1, "[[", 2)
a3<-rep("",20)
a3[16:19]<-c("_Max","_Mean","_Min","_Range")


a4<-paste(unlist(a2),a3,sep="")

names(bios2)<-a4

########## Soil from Soilgrids.com average of each variable for the first 2 meters 

igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
AOI_igh<-st_transform(AOI,igh)


url = "/vsicurl/https://files.isric.org/soilgrids/latest/data/" # Path to the webDAV data.


### Define and download variable of interest

voi = c("bdod","sand","cfvo","clay","phh2o","silt","soc") # variable of interest
depth = c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm","100-200cm")
conv.factor<-c(100, 10,10,10,10,10,10)
quantile = "mean"
(variable = paste(url, voi, sep=""))
(layer = paste(variable,depth,quantile, sep="_")) # layer of interest
(vrt_layer = paste(layer, '.vrt', sep=""))

for (i in 2: length(voi)){
  
  vardepths<-NULL
  
  for (d in 1:length(depth)){
    
    
    datos = paste("/vsicurl/https://files.isric.org/soilgrids/latest/data",voi[i], paste(voi[i],depth[d],"mean.vrt",sep="_"), sep="/")
    paste("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Soil/",paste(voi[i],depth[d],"mean.tif",sep="_"),sep="/")
    
    s1<-terra::rast(datos)
    s2<- terra::crop(s1,ext(vect(AOI_igh)))
    vardepths<-c(vardepths,s2)
    print(paste(voi[i],depth[d]))
  }
  
  vardepths<-rast(vardepths)
  final<-mean(vardepths)
  final1<-terra::project( final,"ESRI:102003", res=2000)
  final2<-terra::focal(
    x = final1,
    w = 7, # sets a 5 cell window around NA pixel
    fun = "mean",
    na.policy = "only", # only interpolate NA values
    na.rm = T,
    overwrite = TRUE
  )
  
  writeRaster(final2,paste("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Soil/",paste(voi[i],"0-2_mean.tif",sep="_"),sep="/"))
  }
  

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Soil/")
soil<-list.files("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Soil/")

soil <- lapply(soil, raster)
soil<-stack(soil)
soil1<-projectRaster(soil,crs = st_crs(AOI)$proj4string)

soil2<-crop(soil1,AOI)
soil2<-mask(soil2,AOI)
soil3<-resample(soil2,bios2[[1]])


test<-stack(soil3,bios2)
####### Sampling bias

plot(st_geometry(AOI[,1]))
plot(st_geometry(bkg1[,1]), add=T)



rr<-owin(xrange=c(-127, -50), yrange=c(8  , 60))
data_rr <- ppp(bkg$decimalLongitude, bkg$decimalLatitude, window = rr)

k04 <- density(data_rr, sigma = 2)

k<-raster(k04)
plot(k)
plot(st_geometry(AOI[,1]),add=T)
k1<-(k+1)
plot(k1)



bias<-crop(k1,AOI)
bias<-mask(bias,AOI)
bias1<-resample(bias,bios2[[1]])
plot(bias1)
names(bias1)<-"Sampling_Bias"

variables<-stack(test,bias1)

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_BIOS")
writeRaster(variables, filename=names(variables), bylayer=TRUE,format="ascii")


vars<-list.files("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_BIOS", pattern=".asc")

vars <- lapply(vars[-c(1,2)], raster)

vars<-stack(vars)
crs(vars)


vars1<-stack(vars,bios2)
vars2<-mask(vars1,AOI)


writeRaster(vars1, filename=names(vars1), bylayer=TRUE,format="ascii")


vars<-list.files("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_BIOS", pattern=".asc")
bios2
for (i in 3:24){
  x<-raster(vars[i])
  x<-projectRaster(x,bios2)
  x<-mask(bios2,AOI)
  
  writeRaster(x, filename=paste("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_BIOS2",names(x),sep="/"),format="ascii")
  
  
}


##################  PCA

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_BIOS2")

vars<-list.files("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_BIOS2")

vars <- lapply(vars[-19], raster)

vars<-stack(vars)

soil.pca <- raster.pca(vars[[c(16:19,21:23)]], 5)
#soil<-raster_pca(vars[[c(16:19,21:23)]])
summary(soil.pca$pca.object)
plot(soil.pca$rasters[[2]])


plot.PCA<-function(pca,pc=c(1,2)){
  par(pty = "s",
      cex.main = 1.2,
      cex.lab = 1,
      font.main = 2,
      font.lab = 2,
      family = "sans",
      col.main = "gray10",
      col.lab = "gray10",
      fg = "gray10",
      las = 1)
  plot.new()
  point.sample<-pca$x[sample(1:length(pca$x[,1]), 10000, replace=F),pc]
  plot.window(xlim = quantile(point.sample[,1],c(0.25,0.75)), 
              ylim =  quantile(point.sample[,2],c(0.25,0.75)), 
              asp = 1)
  
  axis(side = 1, 
       at = round(seq (from=quantile(point.sample[,1],c(0.25)), to=quantile(point.sample[,1],c(0.75)), length.out=3),digits=2),
       labels = TRUE)
  axis(side = 2, 
       at = round(seq(from=quantile(point.sample[,2],c(0.25)), to=quantile(point.sample[,2],c(0.75)), length.out=3),digits = 2),
       labels = TRUE)
  
  title(xlab = paste("PC 1 (", 
                     round(summary(pca)$importance[pc[1]] *100, 
                           digits = 1),
                     "%)", 
                     sep = "") , 
        ylab = paste("PC 2 (", 
                     round(summary(pca)$importance[pc[2]]*100, 
                           digits = 1),
                     "%)", 
                     sep = ""), 
        line = 3,
        adj = 0.5)
  
  
  points(x = pca$x[sample(1:length(pca$x[,1]), 10000, replace=F),pc],
         cex = 0.5, pch=1, col="#99999988")
  
  arrows(x0 = 0, x1 = pca$rotation[,pc[1]], 
         y0 = 0, y1 = pca$rotation[,pc[2]], 
         col = "navy", 
         length = 0.08, 
         lwd = 2,
         angle = 30)
  
  
  text(x = pca$rotation[,pc[1]], y = pca$rotation[,pc[2]], 
       labels = row.names(pca$rotation), 
       cex = 1.2,
       font = 2,
       col = "gray10", 
       pos = c(4, 3, 2, 1, 3, 1))
  
  
}

plot.PCA(soil.pca$pca.object)

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_PCA_BIOS/")
writeRaster(soil.pca$rasters, filename=paste("Soil",names(soil.pca$rasters),sep="_"), bylayer=TRUE,format="ascii")


clim.pca <- raster.pca(vars[[1:15]], 5)
plot(clim.pca$rasters[[2]])
plot.PCA(clim.pca$pca.object)
summary(clim.pca$pca.object)
writeRaster(clim.pca$rasters, filename=paste("Climate",names(clim.pca$rasters),sep="_"), bylayer=TRUE,format="ascii")


###### Western Population

AOI_W<-read_sf("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/AOI_WesternPop.shp")
plot(AOI_W)

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_PCA_BIOS/")

vars<-list.files("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/NA_Current_PCA_BIOS/")

vars <- lapply(vars[-6], raster)

vars<-stack(vars)
plot(vars[[2]])

vars1<-mask(vars,AOI_W)
plot(vars1[[2]])

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Western_Current_PCA/")
writeRaster(vars1, filename=names(vars1), bylayer=TRUE,format="ascii")


### Estern population
AOI_E<-read_sf("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/AOI_EasternPop.shp")
plot(AOI_E)

vars2<-mask(vars,AOI_E)
plot(vars2[[2]])

setwd("C:/Users/seboc/Box/Argasid_tick_project/Maxent_models/Eastern_Current_PCA/")
writeRaster(vars2, filename=names(vars2), bylayer=TRUE,format="ascii")
