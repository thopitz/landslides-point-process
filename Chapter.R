## Due to confidentiality reasons, no real data of landslide occurrences are provided.
## The provided data set has been simulated based on one of the models that showed a good fit to the real data.
## For simplicity, we use only a single working directory for files related to data and output.
## Required R packages are "INLA" (http://www.r-inla.org/download), "spdep", and "fields" (the latter only if maps are plotted).
## Running this script requires files "datamatrix.txt" "SU.shp", "SU.dbf" and "SU.shx" (the latter three only for the spatial effect model), and "Land_Use.txt" (only for plotting Land use categories).

#######################################
############ Preprocessing ############
#######################################

## Set your working directory
setwd("~/research/implementations/landslides/bookchapter/")

## Generate adjacency graph of slope units in file file adjgraph.txt
library(maptools)
SU=readShapeSpatial("SU.shp")

library(spdep)
nb2INLA("adjgraph.txt",poly2nb(SU,queen=F,row.names=SU$ID))

## Load data into data frame: here given in matrix format (one line per pixel)
dataDF=read.delim("datamatrix.txt",header=TRUE,sep=" ")
## Names(dataDF) gives the following variable names:
## Number of landslide events observed at each pixel (integer): status 
## "X" "Y" "status" "Dist2Fault" "NDVI" "RSP" "Slope" "SPI" "TWI" "idsu" "Aspect" "Landform" "Landuse" "Plan_Cur" "Prof_Cur" "Lithology" "Elevation"  
## Description of variables:
## Coordinates in m: X,Y
## Numerical covariates: Dist2Fault,NDVI,RSP, Slope,SPI,TWI,Plan_Cur,Prof_Cur,Elevation
## Categorical covariates: Landform,Landuse,Lithology, Aspect (16 classes)
## Slope unit identifier: idsu
dataDF[dataDF==-9999 | dataDF==-99999]=NA ## Here we set to NA the missing values generated in GIS
data_scaled=dataDF
vars2scale=c("Elevation","Dist2Fault","NDVI","Plan_Cur","Prof_Cur","Slope","SPI","TWI")
data_scaled[,vars2scale]=apply(data_scaled[,vars2scale],2,scale) ## Here we substract the mean and divide by the standard deviation all the continuous covariates

## The subsequent block is used to set the parameters for the model
y.count=data_scaled$count
n.pixels=nrow(data_scaled)
area.pixel=rep(15^2,n.pixels)
offset=area.pixel
n.landslides=sum(y.count)
avg.global=n.landslides/(n.pixels*area.pixel)

## Dataset creation 
covar.inla=data_scaled[,c("Aspect", "Elevation","Dist2Fault","Lithology","Landuse","Landform","NDVI","Plan_Cur","Prof_Cur","Slope","SPI","TWI","SU.ID")]
covar.inla=cbind(intercept=1,covar.inla)
doReplace=covar.inla$Lithology %in% c(8,14,15,16,17,18,20,21,22,23,24,29)
covar.inla$Lithology[doReplace]=0
covar.inla$Aspect=cut(covar.inla$Aspect,labels=F,breaks=seq(0,360,length=17))

####################################################
############ Cox Point Process in INLA  ############
####################################################

library(INLA)
form=y~ -1+intercept+Elevation+Dist2Fault+NDVI+Plan_Cur+Prof_Cur+Slope+SPI+TWI+
     f(Aspect,model="rw1",cyclic=T,constr=T,hyper=list(theta=list(initial=log(5^2),fixed=T)))+
     f(Landuse,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=T)),constr=T)+
     f(Landform,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=T)),constr=T)+
     f(Lithology,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=T)),constr=T)+
     f(SU.ID,model="besag",graph="adjgraph.txt",hyper=list(theta=list(initial=log(1),fixed=F,prior="loggamma",param=c(0.25,0.25))))

stack=inla.stack(data=list(y=y.count,e=offset),A=list(1),effects=list(covar.inla))

fit=inla(form,family="poisson",data=inla.stack.data(stack),
         control.fixed=list(prec=2,prec.intercept=1,mean.intercept=log(avg.global)),
         E=inla.stack.data(stack)$e,num.threads=2)

#########################################
############ Postprocessing  ############
#########################################

## Explore the fitted model ####
## Show summary of fitted model
summary(fit)
     
## If we are interested in looking at significant covariates only, we can use the following function

extractSignificantEffects=function(results.df){
  print(results.df[sign(results.df$`0.025quant`)==sign(results.df$`0.975quant`),])
}

## Show only the significant covariate effects (at 95% level)
extractSignificantEffects(fit$summary.fixed)
extractSignificantEffects(fit$summary.random$Lithology)
extractSignificantEffects(fit$summary.random$Landform)
extractSignificantEffects(fit$summary.random$Landuse)

###########################################
############ Plotting Results  ############
###########################################

## Examples in Line (Aspect) and Point (Land Use) style. 

## Plotting fitted nonlinear effect for Aspect with 95% pointwise confidence envelopes
ylim=range(fit$summary.random$Aspect[,c("mean","0.025quant","0.975quant")])
xvals=((-1):16+0.5)/16*360
yidx=c(16,1:16,1)
plot(xvals, fit$summary.random$Aspect$mean[yidx],type="l",lwd=3,xlab="Angle",ylab="Linear predictor",ylim=ylim,xaxp=c(0,360,9))
lines(xvals, fit$summary.random$Aspect$`0.025quant`[yidx],lwd=3,col="blue")
lines(xvals, fit$summary.random$Aspect$`0.975quant`[yidx],lwd=3,col="blue")
lines(c(0,360),c(0,0),lty=2,lwd=3,col="gray50")

## Plot Landuse results 
lunames=readLines("Land_Use.txt")[-1]
lunames=unlist(strsplit(lunames, "\t"))
lunames=lunames[2*(1:(length(lunames)/2))]
ylim=range(fit$summary.random$Landuse[,c("mean","0.025quant","0.975quant")])
plot(1:13, fit$summary.random$Landuse$mean,lwd=2,xlab="Land-use category",ylab="Linear predictor",ylim=ylim+c(-.12,.0),pch=19,xaxt="n",cex.lab=1)
axis(1,at=as.numeric(rownames(fit$summary.random$Landuse)),labels=rownames(fit$summary.random$Landuse))
points(1:13, fit$summary.random$Landuse$`0.025quant`,lwd=2,col="blue",pch=19)
points(1:13, fit$summary.random$Landuse$`0.975quant`,lwd=2,col="blue",pch=19)
segments(1:13,fit$summary.random$Landuse$`0.025quant`,1:13,fit$summary.random$Landuse$`0.975quant`,lwd=2)
lines(c(0,20),c(0,0),lty=2,lwd=3,col="gray50")
text(1:13-.3,y=fit$summary.random$Landuse$mean,lunames,srt=90,cex=1)

## Plot map of fitted intensities and spatial effect
## Fitted intensity (space unit is m^2)
library(fields)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=fit$summary.fitted.values$mean,nx=500,ny=500)
points(dataDF$X[dataDF$count>0],dataDF$Y[dataDF$count>0],pch=19,cex=.15)
#fitted intensity on log-scale (space unit is km^2)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=log(fit$summary.fitted.values$mean),nx=500,ny=500)
points(dataDF$X[dataDF$count>0],dataDF$Y[dataDF$count>0],pch=19,cex=.15)
vals=fit$summary.random$idsu[dataDF$ID.SU,]
quilt.plot(x=dataDF$X,y=dataDF$X,z=vals$mean,nx=500,ny=500)
