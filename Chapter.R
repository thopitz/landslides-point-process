## This code file accompanies the manuscript "Numerical recipes for landslide spatial prediction by using R-INLA: a step-by-step tutorial" by L. Lombardo, T. Opitz, R. Huser, referred to as LOH in the following.
## Due to confidentiality reasons, NO REAL DATA of landslide occurrences are provided (but covariate data are real).
## The provided data set of pixel-based landslide counts has been simulated using the model fitted in the following code, which showed a good fit to the real data.
## For simplicity, we use only a single working directory for files related to data and output.
## Required R packages are "INLA" (http://www.r-inla.org/download), "spdep", "fields" (only if maps are plotted) and pROC (only if ROC and AUC are calculated).
## Running this script requires files "datamatrix.txt" "SU.shp", "SU.dbf" and "SU.shx" (the latter three only for the spatial effect model), and "Land_Use.txt" (only for plotting Land use categories).

#######################################
############ Preprocessing ############
#######################################

##see Section 3.1 in LOH

## Set your working directory
setwd("~/Set.Your.Local.Directory/")
## Generate adjacency graph of slope units into file adjgraph.txt
library(maptools)
SU=readShapeSpatial("SU.shp")
library(spdep)
nb2INLA("adjgraph.txt",poly2nb(SU,queen=F,row.names=SU$ID))

## Load data into data frame: here given in matrix format (one line per pixel)
dataDF=read.delim("datamatrix.txt",header=TRUE,sep=" ")
## names(dataDF) gives the following variable names:
## "X" "Y" "count" "Aspect" "Dist2Fault" "Elevation" "Landform" "Landuse" "Lithology" "NDVI" "Plan_Cur" "Prof_Cur" "Slope" "SPI" "SU.ID" "TWI"
## Description of variables:
## Number of landslide events observed at each pixel (integer): count 
## Coordinates in meter: X,Y
## Numerical covariates: Aspect,Dist2Fault,Elevation,NDVI,Plan_Cur,Prof_Cur,Slope,SPI,TWI
## Categorical covariates: Landform,Landuse,Lithology
## Slope unit identifier: SU.ID
dataDF[dataDF==-9999 | dataDF==-99999]=NA ## Here we set to NA the missing values generated in GIS
dataDF=na.omit(dataDF) #remove pixels with NA values (at the boundary of the study region)
data_scaled=dataDF #Create a copy of the original data frame
vars2scale=c("Elevation","Dist2Fault","NDVI","Plan_Cur","Prof_Cur","Slope","SPI","TWI")
data_scaled[,vars2scale]=apply(data_scaled[,vars2scale],2,scale) ## Here we scale all continuous covariates (i.e., substract the mean and divide by the standard deviation) 

## The subsequent block defines certain parameters for the model
y.count=data_scaled$count #vector of counts
n.pixels=nrow(data_scaled) #number of pixels in the data set
area.pixel=rep(15^2,n.pixels) #area of pixels (in m)
offset=area.pixel
n.landslides=sum(y.count) #overall number of landslides
avg.global=n.landslides/(n.pixels*area.pixel) #average number of landslides per square meter

## Create dataset of preprocessed covariates
covar.inla=data_scaled[,c("Aspect", "Elevation","Dist2Fault","Lithology","Landuse","Landform","NDVI","Plan_Cur","Prof_Cur","Slope","SPI","TWI","SU.ID")]
covar.inla=cbind(intercept=1,covar.inla)
doReplace=covar.inla$Lithology %in% c(8,14,15,16,17,18,20,21,22,23,24,29)
covar.inla$Lithology[doReplace]=0
covar.inla$Aspect=cut(covar.inla$Aspect,labels=F,breaks=seq(0,360,length=17))

## Several examples of how to plot covariates using R 
library(fields)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=dataDF$Elevation,nx=500,ny=500)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=dataDF$Slope,nx=500,ny=500)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=dataDF$Landform,nx=500,ny=500,breaks=-.5+0:10,nlevel=10)

####################################################
############ Cox Point Process in INLA  ############
####################################################

##see Section 3.2 in LOH

library(INLA)
#Define the model formula, including the specification of prior distributions
form=y~ -1+intercept+Elevation+Dist2Fault+NDVI+Plan_Cur+Prof_Cur+Slope+SPI+TWI+
     f(Aspect,model="rw1",cyclic=T,constr=T,hyper=list(theta=list(initial=log(5^2),fixed=T)))+
     f(Landuse,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=T)),constr=T)+
     f(Landform,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=T)),constr=T)+
     f(Lithology,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=T)),constr=T)+
     f(SU.ID,model="besag",graph="adjgraph.txt",hyper=list(theta=list(initial=log(1),fixed=F,prior="loggamma",param=c(0.25,0.25))))

stack=inla.stack(data=list(y=y.count,e=offset),A=list(1),effects=list(covar.inla))

fit=inla(form,family="poisson",data=inla.stack.data(stack),
         control.fixed=list(prec=2,prec.intercept=1,mean.intercept=log(avg.global)),
         E=inla.stack.data(stack)$e,num.threads=2,
         control.predictor=list(compute=TRUE))

#########################################
############ Postprocessing  ############
#########################################

##see Section 4.1 in LOH

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

## Details on the hyperparameter
fit$summary.hyperpar

## Information on the latent spatial effect, its 97.5% an 2.5% posterior percentiles 
fit$summary.random$SU.ID$mean
fit$summary.random$SU.ID$"0.025quant"
fit$summary.random$SU.ID$"0.975quant"

###########################################
############ Plotting Results  ############
###########################################

##see Section 4.1 in LOH

## Examples in Line (Aspect) and Point (Land Use) style. 

## Plot fitted nonlinear effect for Aspect with 95% pointwise confidence envelopes
aspect=fit$summary.random$Aspect
ylim=range(aspect[,c("mean","0.025quant","0.975quant")])
xvals=(c(-1:16)+0.5)/16*360; yidx=c(16,1:16,1)
plot(xvals, aspect$mean[yidx],type="l",lwd=3,xlab="Aspect [Deg]",ylab="Linear predictor",ylim=ylim,xaxp=c(0,360,9))
lines(xvals, aspect$"0.025quant"[yidx],lwd=3,col="blue")
lines(xvals, aspect$"0.975quant"[yidx],lwd=3,col="blue")
abline(h=0,lty=2,lwd=3,col="gray50")

## Plot Landuse results 
LUnames=readLines("Land_Use.txt")[-1]
LUnames=unlist(strsplit(LUnames, "\t"))
LUnames=LUnames[2*(1:(length(LUnames)/2))]
landuse=fit$summary.random$Landuse
ylim=range(landuse[,c("mean","0.025quant","0.975quant")])
plot(1:13,landuse$mean,xlab="CORINE Land Use",ylab="Linear predictor",lwd=2,ylim=ylim+c(0,0.1),pch=19,xaxt="n")
axis(1,at=1:13,labels=1:13)
points(1:13,landuse$"0.025quant",lwd=2,col="blue",pch=19)
points(1:13,landuse$"0.975quant",lwd=2,col="blue",pch=19)
segments(1:13,landuse$"0.025quant",1:13,landuse$"0.975quant",lwd=2)
abline(h=0,lty=2,lwd=3,col="gray50")
text(1:13-0.3,y=landuse$mean,LUnames,srt=90)

## Plot map of fitted intensities and spatial effect
## Fitted pixel intensity (space unit is m^2)
library(fields)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=fit$summary.fitted.values$mean,nx=500,ny=500)
points(dataDF$X[dataDF$count>0],dataDF$Y[dataDF$count>0],pch=19,cex=.15)
## Fitted pixel intensity on log-scale (space unit is m^2)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=log(fit$summary.fitted.values$mean),nx=500,ny=500)
points(dataDF$X[dataDF$count>0],dataDF$Y[dataDF$count>0],pch=19,cex=.15)
## Fitted latent spatial effect (space unit is m^2)
vals=fit$summary.random$SU.ID[dataDF$SU.ID,]
quilt.plot(x=dataDF$X,y=dataDF$Y,z=vals$mean,nx=500,ny=500)
## Size of 95% credible interval of latent spatial effect (space unit is m^2)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=vals$`0.975quant`-vals$`0.025quant`,nx=500,ny=500)
## Significance of  latent spatial effect (space unit is m^2)
slow=sign(vals$`0.025quant`)
sup=sign(vals$`0.975quant`)
vals.signif=ifelse(slow==sup,ifelse(slow<0,-1,1),NA)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=vals.signif,nx=500,ny=500,nlevel=2,add.legend=F)
## Fitted SU intensity
pixel.intensity=fit$summary.fitted.values$mean*area.pixel
SU.intensity=aggregate(pixel.intensity,by=list(SU.ID=covar.inla$SU.ID),FUN=sum)$x
quilt.plot(x=dataDF$X,y=dataDF$Y,z=SU.intensity[dataDF$SU.ID],nx=500,ny=500)

## Examples of Spatial plots for fitted covariate effects
library(fields)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=fit$summary.random$Aspect$mean[covar.inla$Aspect],nx=500,ny=500)
quilt.plot(x=dataDF$X,y=dataDF$Y,z=fit$summary.random$Landform$mean[covar.inla$Landform+1],nx=500,ny=500)


##################################################
##### Model validation and cross-validation  #####
##################################################

##see Sections 4.2 and 4.3 in LOH

## Notice that this code section is very time- and memory-intensive since it requires K=4 full runs of the INLA fit.

## Randomly partition slope units into K=4 subsections
K=4
n.SU=max(covar.inla$SU.ID)
partSU=sample(as.integer(cut(1:n.SU,breaks=K)))
## Create partition at the pixel level
partition=partSU[covar.inla$SU.ID]
fits=list() # List to store all cross-validation fits
pred.intensity=c() # Vector to store all predicted pixel intensities
for(i in 1:K){
  sec.i=which(partition==i) # Pixel indices for the i-th subsection
  y.count.i=y.count # Clone response vector
  y.count.i[sec.i]=NA # Hide i-th subsection for fitting
  stack.i=inla.stack(data=list(y=y.count.i,e=offset),A=list(1),effects=list(covar.inla)) # Create new data stack
  fit.i=inla(form,family="poisson",data=inla.stack.data(stack.i),
           control.fixed=list(prec=2,prec.intercept=1,mean.intercept=log(avg.global)),
           E=inla.stack.data(stack.i)$e,num.threads=2,
           control.predictor=list(compute=TRUE,link=1)) # Perform the i-th INLA fit
  fits[[i]]=fit.i # Store the i-th cross-validation fit
  pred.intensity[sec.i]=fit.i$summary.fitted.values$mean[sec.i]*area.pixel 
}

## Prepare observed and fitted values at pixel and slope unit scale (within-sample and cross-validation)
pixel.intensity=fit$summary.fitted.values$mean*area.pixel
pixel.intensity.cv=pred.intensity
SU.intensity=aggregate(pixel.intensity,by=list(SU.ID=covar.inla$SU.ID),FUN=sum)$x
SU.intensity.cv=aggregate(pixel.intensity.cv,by=list(SU.ID=covar.inla$SU.ID),FUN=sum)$x # For cross-validation
SU.count=aggregate(y.count,by=list(SU.ID=covar.inla$SU.ID),FUN=sum)$x

## Plot observed vs. fitted intensities at the slope unit level
plot(SU.count,SU.intensity, pch=19, xlim = c(0,60), ylim = c(0,60), cex=0.6, xlab = "Observed Landslide Count", ylab = "Estimated Intensity", main = "Slope Units", font.main = 1)
points(SU.count, SU.intensity.cv, pch=19, cex=0.6)
abline(0,1, col = "red")
lines(seq(0,70,length=10000),qpois(0.025,lambda=seq(0,70,length=10000)),col="lightgrey")
lines(seq(0,70,length=10000),qpois(0.975,lambda=seq(0,70,length=10000)),col="lightgrey")
legend("bottomright",legend=c("Fit","CV"),col=c("cyan3","black"),pch=19)


## Calculate ROC curves and AUC values at the pixel scale
library(pROC)
y.presence.absence=y.count
y.presence.absence[y.count>0]=1
pixel.probability=1-exp(-pixel.intensity) 
pixel.probability.cv=1-exp(-pixel.intensity.cv) 
SU.presence.absence=SU.count
SU.presence.absence[SU.count>0]=1
SU.probability=1-exp(-SU.intensity)
SU.probability.cv=1-exp(-SU.intensity.cv)

pixel.ROC=roc(y.presence.absence~pixel.probability)
pixel.ROC.cv=roc(y.presence.absence~pixel.probability.cv)
pixel.AUC=as.numeric(pixel.ROC$auc)
pixel.AUC.cv=as.numeric(pixel.ROC.cv$auc)

SU.ROC=roc(SU.presence.absence~SU.probability)
SU.ROC.cv=roc(SU.presence.absence~SU.probability.cv)
SU.AUC=as.numeric(SU.ROC$auc)
SU.AUC.cv=as.numeric(SU.ROC.cv$auc)


## Plot pixel ROC 
plot(1-pixel.ROC$specificities,pixel.ROC$sensitivities,xlim=c(0,1),ylim=c(0,1),col="cyan3",xlab="1-Specificity",ylab="Sensitivity",main="ROC curve")
lines(1-pixel.ROC.cv$specificities,pixel.ROC.cv$sensitivities)
legend("bottomright",legend=c("Fit","CV"),col=c("cyan3","black"),lty=c(1,1))

## Plot SU ROC
plot(1-SU.ROC$specificities,SU.ROC$sensitivities,xlim=c(0,1),ylim=c(0,1),col="cyan3",xlab="1-Specificity",ylab="Sensitivity",main="ROC curve")
lines(1-SU.ROC.cv$specificities,SU.ROC.cv$sensitivities)
legend("bottomright",legend=c("Fit","CV"),col=c("cyan3","black"),lty=c(1,1))




