# Load the Libraries
install.packages('gganimate')
##############################
library(ggplot2)
library(dplyr);library(viridis)
library(ggExtra)
library(magrittr)
library(gstat)
library(phylin)
library(sp) 
library(tidyverse)
library(sf)
library(ape)
library(gstat)
library(spatstat)
library(Metrics)
library(spdep)
library("plot3Drgl")
library(spatial)
library(gganimate)

#Step 1:Explanatory Data Analysis

#Reading the data in R, Viewing it and summarizing it.
Data=read.csv("Temp.csv")
View(Data)
summary(Data)


mycol <- colorRampPalette(c("#eeaf61", "#fb9062", "#ee5d6c", "#ce4993", "#6a0d83"))

mycol2 <- colorRampPalette(c("#a4e7da", "#fff2cc" , "#ffad7d", "#ff6b6b", "#ff2a2a"))

#Box Plot
ggplot(data=Data,aes(y=Average,x=""))+
  geom_boxplot(fill = "#ff6b6b",alpha=0.3,outlier.colour ="#ee5d6c",colour="#fb9062")+ 
  ggtitle("Temperature Boxplot") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#Another box plot function
#boxplot(Data$Average)

# Create a histogram of average temperature
ggplot(data=Data, aes(x =Average)) +
  geom_histogram(binwidth = 1, fill = "#ee5d6c") +
  ggtitle("Average Temperature Histogram") +
  xlab("Average Temperature") +
  ylab("Frequency")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#Another Histogram Function
#hist(Data$Average)


# Global Moran I
temp.dists<-as.matrix(dist(cbind(Data$LONG,Data$LAT)))
temp.dist.inv<-1/temp.dists
diag(temp.dist.inv)<-0

l.moransI=Moran.I(Data$Average,temp.dist.inv)
l.moransI


#Moran
#rule of thumb for selecting the number of neighbors in KNN analysis
#is to use the square root of the number of observations in your data set
library(lctools)
l.moran<-l.moransI(temp.dists,20,Data$Average, WType='Bi-square', scatter.plot = TRUE, family = "adaptive")
l.moran


# lisa plot 
library(ncf)
lisa <- lisa(Data$LONG, Data$LAT, Data$Average, neigh = 20, latlon=TRUE)
lisa
plot(lisa)

#Ask Dr sara Elsob7

#bubble plot (1)
symbols(Data$LONG, Data$LAT, circles=Data$Average)

#bubble plot (2)
plot(Data$LONG, Data$LAT,colour=Data$Average)

#bubble plot(3)
# Create a scatter plot with points colored by "average"
ggplot(Data, aes(x = LONG, y = LAT, color = Average)) +
  # Add points with size proportional to a third variable
  geom_point(aes(size = Average)) +
  # Add a legend for the color scale
  scale_color_gradient(low = "blue", high = "red")

#bubble plot(4)
ggplot(Data, aes(x = LONG, y = LAT, color = Average)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red", 
                       na.value = "gray", 
                       name = "Average")

#bubble plot(5) A7laaa wa7edd
ggplot(Data, aes(x = LONG, y = LAT, color = Average)) +
  geom_point(aes(size = Average)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


#3D plots
scatter3D(Data$LONG,Data$LAT,Data$Average, zcol=Data$Average)
scatter3D(Data$LONG,Data$LAT,Data$Average, zcol=Data$Average,ticktype="detailed")
scatter3D(Data$LONG,Data$LAT,Data$Average, zcol=Data$Average,pty="g",ticktype="detailed")
#############################################################################################
###Spatial Interpolation

##Fitting a trend surface model by least squares
library(spatial)
x<-Data$LONG
y<-Data$LAT
z<-Data$Average

#3rd degree

  fit.sfc3 <- surf.ls(3,x,y,z)
summary(fit.sfc3)
fit.sfc3$beta

#2nd degree

fit.sfc2 <- surf.ls(2,x,y,z)
summary(fit.sfc2)
fit.sfc2$beta

#4th  degree

fit.sfc2 <- surf.ls(4,x,y,z)
summary(fit.sfc2)
fit.sfc2$beta

#5th  degree

fit.sfc2 <- surf.ls(5,x,y,z)
summary(fit.sfc2)
fit.sfc2$beta


#Compare between them using AIC, choose the smallest AIC evaluate trend surface over a grid
trsurf3 <- trmat(fit.sfc2, min(x), max(x), min(y), max(y), 50)

trsurf3

scatter3D(x,y,z,surf = list(x = trsurf3$x, y = trsurf3$y, z = trsurf3$z,
                            NAcol = "grey", shade = 0.1))


library(plot3Drgl)

contour(trsurf3)
##################################################################
## IDW

###IDW
library(phylin)
library(sp) 

#Extracting the columns we need from our data (Long, Lat and Average column)
#We want to put the 3 columns in a data frame

#Way 1

coords <- Data[ , c("LONG", "LAT")]

data1  <- Data[ , c("Average")]

#Way 2


coords <- Data[ , c("LONG", "LAT")]   # coordinates
data   <- as.data.frame(Data[ , c("Average")])          # data


#This Line of the code used for projection
crs    <- CRS("+init=epsg:28992") # proj4string of coords

# make the SpatialPointsDataFrame object
spdf <- SpatialPointsDataFrame(coords      = coords,
                               data        = data, 
                               proj4string = crs)


#Making Interpolated maps

Long <- seq(from=-124.23,to=-114.16,by=0.1)
Lat <- seq(from=32.55,to=41.98,by=0.1)
grid<-cbind(rep(Long,length(Lat)), rep(Lat,each=length(Long)))

idw<- phylin::idw(data1, coords, grid,p=3.991)


grid.image(idw, grid, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2")
points(coords, cex=data1/100)

grid.image(idw, grid, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = rainbow)
points(coords, cex=data1/100)

grid.image(idw, grid, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab=paste ("Genetic distance to sample s2"), colFUN = terrain.colors)
points(coords, cex=data1/100)

mycolors <- colorRampPalette(c("#6a0dad", "#9e36b2", "#d163be", "#f89ece", "#b1dafb"))
grid.image(idw, grid, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = mycol)


mycolors1 <- colorRampPalette(c("#f5a9b8", "#f9f196", "#86e9b4", "#81d7f7", "#d996f8"))
grid.image(idw, grid, main='IDW interpolation', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = mycolors1)

#############################################################################
#IDW OPTIMUM 
library(spatstat)

obs_window <- owin(xrange = c(-124.23,-114.16), yrange =c(32.55,41.98) )

ppp_av<-ppp(Data$LONG,Data$LAT,
            marks=Data$Average,window=obs_window)

idw_av <- idw(ppp_av, power=3.991, at="pixels")

plot(idw_av,
     col=heat.colors(20), 
     main="Interpolated spatial variation in risk of malaria based on IDW method \n (Power = 0.05)") 

idw_points <- idw(ppp_av, power = 0.05, at = "points")

Metrics::mse(ppp_av$marks, idw_points)

powers <- seq(0.001, 10, 0.01)
mse_result <- NULL
for(power in powers){
  CV_idw <- idw(ppp_av, power=power, at="points")
  mse_result <- c(mse_result,
                  Metrics::mse(ppp_av$marks,CV_idw))
}
optimal_power <- powers[which.min(mse_result)]
optimal_power

plot(powers, mse_result)

#################################################################################
######Kriging
# Calculate the empirical variogram of the 'Average' variable 
# The 'Average~1' formula specifies the variable to be analyzed and '1' indicates that there are no covariates
# The 'boundaries' argument sets the distance lags for the variogram calculation
evgm <- variogram(Data[, c("Average")]~1,spdf)
plot(evgm)

# Fit a variogram model to the experimental variogram
# The 'vgm' function specifies the type of variogram model to be used
# Here, a Matérn model with a sill of 5, a range of 0.7 and a nugget of 1 is used
fvgm1 <- fit.variogram(evgm,vgm(5,"Mat",0.7,1))
fvgm1

fvgm2 <- fit.variogram(evgm,vgm("Sph"))
fvgm2
# Plot the experimental variogram with the fitted variogram model
# The 'model' argument specifies the fitted variogram model to be plotted
plot(evgm,model=fvgm1)
plot(evgm,model=fvgm2)

#extract the value of the sum of squared errors from the fitted variogram models 
#fvgm1 and fvgm2, respectively. The sum of squared errors is a measure 
#of the difference between the observed values and the predicted values from the model.
#A smaller value of SSE indicates a better fit of the model to the data.
attr(fvgm1, "SSErr")
attr(fvgm2, "SSErr")


#GAUSSIAN
fvgm3 <- fit.variogram(evgm,vgm("Gau"))
fvgm3
plot(evgm,model=fvgm3)
attr(fvgm3, "SSErr")
#EXPONENTIAL
fvgm4 <- fit.variogram(evgm,vgm("Exp"))
fvgm4

summary(fvgm1)
plot(evgm,model=fvgm4)
attr(fvgm4, "SSErr")

# fitting simple kriging on a grid
s.grid <- spsample(spdf, type = "regular", n = 6000)

krig.est <- krige(data[, c("Average")] ~ 1, spdf, newdata = s.grid, model = fvgm1)

spplot(krig.est["var1.pred"])
spplot(krig.est["var1.var"])
#changing range
library(sp)
data(data)
coordinates(temp) <- LONG+LAT
data(krig.est)
gridded(krig.est) <- LONG+LAT
spplot(krig.est["dist"], at = seq(0, 1000, by = 100))
#############

krig.grid <- SpatialPixelsDataFrame(krig.est, krig.est@data)
levs <- c(0, 20, 50, 70, 100, Inf)
var.levs <- c(0, 150, 190, 220, 250, Inf)
#using t map but proplem in variance map
library(tmap)
krig.map.est <- tm_shape(krig.grid) +
  tm_raster(col = 'var1.pred', breaks = levs, title = 'data', palette = 'Reds') +
  tm_layout(legend.bg.color = 'white', legend.frame = TRUE)

krig.map.var <- tm_shape(krig.grid) +
  tm_raster(col = 'var1.var', breaks = var.levs, title = 'Estimate Variance', palette = 'Reds') +
  tm_layout(legend.bg.color = 'white', legend.frame = TRUE)

tmap_arrange(krig.map.est, krig.map.var)


####### cross validation########
##### remember that you can use training data#######



coordinates(ozone)<-~Lon+Lat
proj4string(ozone) <- crs

krig.est <- krige(Av8top ~ 1,  ozone, newdata = s.grid, model = fvgm2)
cv.o <- krige.cv(Av8top ~ 1, ozone, fvgm2, nfold=nrow(ozone))
summary(cv.o)

res <- as.data.frame(cv.o)$residual
sqrt(mean(res^2))

mean(res)

mean(res^2/as.data.frame(cv.o)$var1.var)

##### using grid.image library((phylin)#########
summary(spdf)
library(phylin)
Long <- seq(from=-125,to=-114,by=0.01)
Lat <- seq(from=30,to=42,by=0.01)
grid<-expand.grid(Long,Lat)
gridsp<- SpatialPoints(grid,proj4string = crs)


krig.est <- krige(data[, c("Average")] ~ 1, spdf, newdata = gridsp, model = fvgm1)
kriging<-as.matrix(krig.est$var1.pred)
grid.image(kriging, grid, main='ordinary kriging', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = terrain.colors)


##### interpolated map########

grid.image(kriging, grid, main='ordinary kriging', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2", colFUN = terrain.colors)
points(coords, cex=data1/3)

#######prediction error map##########

variance<-as.matrix(krig.est$var1.var)
grid.image(variance, grid, main='ordinary kriging', xlab='Longitude', 
           ylab='Latitude', sclab="Genetic distance to sample s2")
points(coords, cex=data1/3)
##############

#geographical weighted regression model
library(spgwr)
library(ggplot2)
model1 <- lm(Data$Average ~ Data$LONG + Data$LAT)
summary(model1)
plot(model1)
par(mfrow=c(2,2))
plot(model1)


library(spatial)
x<-Data$LONG
y<-Data$LAT
z<-Data$Average


resids<-residuals(model1)
colours <- c("dark blue", "blue", "red", "dark red")
map.resids <- SpatialPointsDataFrame(data = data.frame(resids), coords = cbind(x,y))
spplot(map.resids, cuts=quantile(resids), col.regions=colours, cex=1) 


#calculate kernel bandwidth
GWRbandwidth <- gwr.sel(Data$Average ~ Data$LAT+Data$LONG, data=Data, coords=cbind(x,y),adapt=T) 

gwr.model = gwr(Data$Average ~ Data$LAT+Data$LONG, data=Data, coords=cbind(x,y), adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 
gwr.model
results<-as.data.frame(gwr.model$SDF)
head(results)

Data$LONG<-results$LONG
Data$LAT<-results$LAT
gwr.point1<-ggplot(Data, aes(x=x,y=y))+geom_point(aes(colour=x))+scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Coefs"))
gwr.point3<-ggplot(Data, aes(x=x,y=y))+geom_point(aes(colour=y))+scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Coefs"))
gwr.point3+geom_path(Data,aes(x, y, group=id), colour="grey")+coord_equal()
plot(gwr.point1)
plot(gwr.point3)
plot(gwr.point3+geom_path(Data,aes(x, y, group=id), colour="grey")+coord_equal())





