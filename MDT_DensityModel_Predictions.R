#===============================================================================================# 

 ## Authors: ER Zylstra, RJ Steidl, N Pope

 ## Project: MDT line-distance sampling - estimating density (spatial variation and trends)

  # Model structure adapted from Chelgren et al. 2011, Herpetol. Conserv. Biol. 6:175-190
  # and Shirk et al. 2014, Diversity and Distributions 20:1186-1199

 ## Relevant files: 

  # 1. MDT_Import_Format_Data.R -- import all the files needed for analysis, format data, and create MDT_Data.Rdata
  # 2. MDT_DensityModel.R -- load MDT_Data.Rdata, identify covariates to include in model, run model in STAN, and view results
  # 3. MDT_DensityModel.stan -- STAN model
  # 4. MDT_Predictions.R -- predict tortoise density throughout the study area

#===============================================================================================#  

#-----------------------------------------------------------------------------------------------# 
# Set working directory, load packages
#-----------------------------------------------------------------------------------------------# 

  #Set working directory

    #setwd()

  #Load packages
    
    library(plyr)
    library(reshape2)
    library(raster)
    library(rgdal)
    library(rstan)

  # rm(list=ls())
  
#-----------------------------------------------------------------------------------------------# 
# Import data
#-----------------------------------------------------------------------------------------------# 
  
  #Import survey and associated covariate data (so covariate rasters can be standardized appropriately)

    load('MDT_Data.Rdata')
    
  #Output from Stan model (and associated covariate data that went into the model)
    
    load('xx.Rdata')  #assuming the Stan fitted object (out) was saved in .Rdata file  
    posterior6000 <- as.matrix(out)
    
  #Covariate data
    
    precip.s <- raster('Covariates/30yr_may_oct_prcp_mm.tif')            #Summer precip (May-Oct, mm), 30yr norms
    names(precip.s) <- 'precips'
    precip.w <- raster('Covariates/30yr_nov_apr_prcp_mm.tif')            #Winter precip (Nov-Apr, mm), 30yr norms
    names(precip.w) <- 'precipw'
    elev <- raster('Covariates/avg_elevation.tif')                       #Elevation (m)
    names(elev) <- 'elev'
    north <- raster('Covariates/avg_northness.tif')                      #Aspect:northness 
    names(north) <- 'north'
    east <- raster('Covariates/avg_eastness.tif')                        #Aspect:eastness 
    names(east) <- 'east'
    rough <- raster('Covariates/average_surface_roughness_snapped.tif')  #Average surface roughness
    names(rough) <- 'rough'
    wash <- raster('Covariates/wash_proportion.tif')                     #Wash density (proportion of each pixel predicted to be "wash")
    names(wash) <- 'wash'
    bedrock <- raster('Covariates/bedrock_depth.tif')                    #Depth to bedrock (cm)
    names(bedrock) <- 'bedrock'
    veg.a <- raster('Covariates/annProx_snapped.tif')                    #Annual plant growth potential (MODIS, compare wet(2005)-dry(2002) years)
    names(veg.a) <- 'vega'
    road.a <- raster('Covariates/dist_road.tif')                         #Distance to any road (m)
    names(road.a) <- 'roada'
    #Note: rough and veg.a are missing values in the far western part of WM (around Mojave CA)
    
  #Shapefile: Recovery units
    rus <- readOGR(dsn='Covariates/Revised Recovery Units',layer='2011RecoveryUnits')
    rus <- spTransform(rus,crs(east))
    
  #Shapefile: State boundaries
    states <- readOGR(dsn='Covariates/US_states_GIS',layer='cb_2017_us_state_500k')
    states <- spTransform(states,crs(east))
    
  #Shapefiles: major roads (keep commented out, except when needed since it's slow)
    # roads.nv <- readOGR(dsn='Covariates/Roads_NV',layer='tl_2021_32_prisecroads')
    # roads.nv <- spTransform(roads.nv,crs(east))
    # roads.ca <- readOGR(dsn='Covariates/Roads_CA',layer='tl_2021_06_prisecroads')
    # roads.ca <- spTransform(roads.ca,crs(east))   
    # roads.az <- readOGR(dsn='Covariates/Roads_AZ',layer='tl_2021_04_prisecroads')
    # roads.az <- spTransform(roads.az,crs(east))  
    # roads.ut <- readOGR(dsn='Covariates/Roads_UT',layer='tl_2021_49_prisecroads')
    # roads.ut <- spTransform(roads.ut,crs(east))

#-----------------------------------------------------------------------------------------------# 
# Format survey data and obtain means/SDs used to standardize covariates associated with each survey
#-----------------------------------------------------------------------------------------------# 

    #Covariates that were log-transformed before standardizing: rough, wash, bedrock, road.a
    
      surveys$rough.log <- log(surveys$rough)
      surveys$wash.log <- log(surveys$wash)
      surveys$bedrock.log <- log(surveys$bedrock)
      surveys$road.a.log <- log(surveys$road.a)
    
    #Getting means, SDs for each covariate in the density model

      precips.mn <- mean(surveys$precip.s)
      precipw.mn <- mean(surveys$precip.w)
      elev.mn <- mean(surveys$elev)
      rough.l.mn <- mean(surveys$rough.log)
      east.mn <- mean(surveys$aspect.e)
      north.mn <- mean(surveys$aspect.n)
      wash.l.mn <- mean(surveys$wash.log)
      bedrock.l.mn <- mean(surveys$bedrock.log)
      vega.mn <- mean(surveys$veg.a)
      roada.l.mn <- mean(surveys$road.a.log)
      precipw.sd <- sd(surveys$precip.w)
      precips.sd <- sd(surveys$precip.s)
      elev.sd <- sd(surveys$elev)
      rough.l.sd <- sd(surveys$rough.log)
      east.sd <- sd(surveys$aspect.e)
      north.sd <- sd(surveys$aspect.n)
      wash.l.sd <- sd(surveys$wash.log)
      bedrock.l.sd <- sd(surveys$bedrock.log)
      vega.sd <- sd(surveys$veg.a)
      roada.l.sd <- sd(surveys$road.a.log)
    
#-----------------------------------------------------------------------------------------------# 
# Transform (where needed) and standardize values in covariate rasters
#-----------------------------------------------------------------------------------------------# 

    precips.z <- (precip.s - precips.mn)/precips.sd
    precipw.z <- (precip.w - precipw.mn)/precipw.sd
    elev.z <- (elev - elev.mn)/elev.sd
    east.z <- (east - east.mn)/east.sd
    north.z <- (north - north.mn)/north.sd
    vega.z <- (veg.a - vega.mn)/vega.sd
    
    #Rough: remove values outside of 1-1.3 range (for now)
      summary(surveys$rough) #survey values range from 1.001 to 1.193
      sum(values(rough)<1 | values(rough)>1.3, na.rm=TRUE) #258 cells have these extreme values
      # plot(rough>1.3)
      rough[rough<1 | rough>1.5] <- NA
      rough.z <- (log(rough) - rough.l.mn)/rough.l.sd
        names(rough.z) <- 'rough.l'
    
    #Bedrock (survey values range from 102 to 48905)
      # plot(bedrock < 102) #probably a good idea to remove cells with values below minimum survey value (some <1)
      bedrock[bedrock<102] <- NA
      bedrock.z <- (log(bedrock) - bedrock.l.mn)/bedrock.l.sd
        names(bedrock.z) <- 'bedrock.l'
    
    #Road (survey values range from 0.13 to 8,650])
      roada.z <- (log(road.a) - roada.l.mn)/roada.l.sd
        names(roada.z) <- 'roada.l'
    
    #Wash
      #Missing some values in the north-central part of the study area
      #Will impute the mean value, like we did before extracting values for survey locations
      wash.missing <- wash
      wash[is.na(wash.missing)] <- mean(values(wash.missing),na.rm=TRUE)  
      #Adjusting minimum value of wash, which was very small, but >0 for all surveys
      wash[wash < min(surveys$wash)] <- min(surveys$wash)
      wash.z <- (log(wash) - wash.l.mn)/wash.l.sd
        names(wash.z) <- 'wash.l'

#-----------------------------------------------------------------------------------------------# 
# Crop covariate layers and divide into smaller units
#-----------------------------------------------------------------------------------------------# 

  #Polygon for each RU
    
    cd <- subset(rus,Unit_Name=='Colorado Desert')
    wm <- subset(rus,Unit_Name=='Western Mojave')
    em <- subset(rus,Unit_Name=='Eastern Mojave')
    nm <- subset(rus,Unit_Name=='Northeastern Mojave')
  
  #Create east.z layer for each recovery unit
    
    east.z.cdA1 <- crop(east.z,cd)
    east.z.cdA <- mask(east.z.cdA1,cd)
    east.z.wmA1 <- crop(east.z,wm)
    east.z.wmA <- mask(east.z.wmA1,wm)
    east.z.emA1 <- crop(east.z,em)
    east.z.emA <- mask(east.z.emA1,em)  
    east.z.nm1 <- crop(east.z,nm)
    east.z.nm <- mask(east.z.nm1,nm)
  
  #Crop Eastern Mojave (exclude areas north of 4060000 and west of 550000)
    
    em.bb <- extent(east.z.emA)
    em.bb[1] <- 550000
    em.bb[4] <- 4060000
    em.poly <- as(em.bb,'SpatialPolygons')
    proj4string(em.poly) <- crs(east.z.emA)
    em.poly <- spTransform(em.poly,proj4string(east.z.emA))
    east.z.em1 <- crop(east.z.emA,em.poly)
    east.z.em <- mask(east.z.em1,em.poly)
    # plot(east.z.em)

  #Crop Western Mojave (exclude areas north of 3950000) and split in half (at easting = 490000)
    
    wmW.bb <- wmE.bb <- extent(east.z.wmA)
    wmW.bb[4] <- wmE.bb[4] <- 3950000
    wmW.bb[2] <- 490000
    wmE.bb[1] <- 490000
    wmW.poly <- as(wmW.bb,'SpatialPolygons')
    wmE.poly <- as(wmE.bb,'SpatialPolygons')
    proj4string(wmW.poly) <- crs(east.z.wmA)
    proj4string(wmE.poly) <- crs(east.z.wmA)
    wmW.poly <- spTransform(wmW.poly,proj4string(east.z.wmA))
    wmE.poly <- spTransform(wmE.poly,proj4string(east.z.wmA))
    east.z.wmW1 <- crop(east.z.wmA,wmW.poly)
    east.z.wmE1 <- crop(east.z.wmA,wmE.poly)
    east.z.wmW <- mask(east.z.wmW1,wmW.poly)
    east.z.wmE <- mask(east.z.wmE1,wmE.poly)

  #Split Colorado Desert in half (northing = 3770000)

    cdN.bb <- cdS.bb <- extent(east.z.cdA)
    cdN.bb[3] <- 3770000
    cdS.bb[4] <- 3770000
    cdN.poly <- as(cdN.bb,'SpatialPolygons')
    cdS.poly <- as(cdS.bb,'SpatialPolygons')
    proj4string(cdN.poly) <- crs(east.z.cdA)
    proj4string(cdS.poly) <- crs(east.z.cdA)
    cdN.poly <- spTransform(cdN.poly,proj4string(east.z.cdA))
    cdS.poly <- spTransform(cdS.poly,proj4string(east.z.cdA))    
    east.z.cdN1 <- crop(east.z.cdA,cdN.poly)
    east.z.cdS1 <- crop(east.z.cdA,cdS.poly)
    east.z.cdN <- mask(east.z.cdN1,cdN.poly)
    east.z.cdS <- mask(east.z.cdS1,cdS.poly)  
    
  #Check (merge and plot): 
    # r6 <- list(east.z.em,east.z.wmW,east.z.wmE,east.z.cdN,east.z.cdS,east.z.nm)
    # east.z6 <- do.call(merge,r6)
    # plot(east.z6)
    # plot(rus,add=TRUE)

  #Crop other rasters to match the east.z layers

    #Colorado Desert - north half
    precips.z.cdN1 <- crop(precips.z,east.z.cdN)
    precipw.z.cdN1 <- crop(precipw.z,east.z.cdN)
    elev.z.cdN1    <- crop(elev.z,east.z.cdN)
    north.z.cdN1   <- crop(north.z,east.z.cdN)
    rough.z.cdN1   <- crop(rough.z,east.z.cdN)
    wash.z.cdN1    <- crop(wash.z,east.z.cdN)
    bedrock.z.cdN1 <- crop(bedrock.z,east.z.cdN)
    vega.z.cdN1    <- crop(vega.z,east.z.cdN)
    roada.z.cdN1   <- crop(roada.z,east.z.cdN)
    precips.z.cdN <- mask(precips.z.cdN1,east.z.cdN)
    precipw.z.cdN <- mask(precipw.z.cdN1,east.z.cdN)
    elev.z.cdN <- mask(elev.z.cdN1,east.z.cdN)
    north.z.cdN <- mask(north.z.cdN1,east.z.cdN)
    rough.z.cdN <- mask(rough.z.cdN1,east.z.cdN)
    wash.z.cdN <- mask(wash.z.cdN1,east.z.cdN)
    bedrock.z.cdN <- mask(bedrock.z.cdN1,east.z.cdN)
    vega.z.cdN <- mask(vega.z.cdN1,east.z.cdN)
    roada.z.cdN <- mask(roada.z.cdN1,east.z.cdN)
    
    #Colorado Desert - south half
    precips.z.cdS1 <- crop(precips.z,east.z.cdS)
    precipw.z.cdS1 <- crop(precipw.z,east.z.cdS)
    elev.z.cdS1    <- crop(elev.z,east.z.cdS)
    north.z.cdS1   <- crop(north.z,east.z.cdS)
    rough.z.cdS1   <- crop(rough.z,east.z.cdS)
    wash.z.cdS1    <- crop(wash.z,east.z.cdS)
    bedrock.z.cdS1 <- crop(bedrock.z,east.z.cdS)
    vega.z.cdS1    <- crop(vega.z,east.z.cdS)
    roada.z.cdS1   <- crop(roada.z,east.z.cdS)
    precips.z.cdS <- mask(precips.z.cdS1,east.z.cdS)
    precipw.z.cdS <- mask(precipw.z.cdS1,east.z.cdS)
    elev.z.cdS <- mask(elev.z.cdS1,east.z.cdS)
    north.z.cdS <- mask(north.z.cdS1,east.z.cdS)
    rough.z.cdS <- mask(rough.z.cdS1,east.z.cdS)
    wash.z.cdS <- mask(wash.z.cdS1,east.z.cdS)
    bedrock.z.cdS <- mask(bedrock.z.cdS1,east.z.cdS)
    vega.z.cdS <- mask(vega.z.cdS1,east.z.cdS)
    roada.z.cdS <- mask(roada.z.cdS1,east.z.cdS)
    
    #Western Mojave - west half
    precips.z.wmW1 <- crop(precips.z,east.z.wmW)
    precipw.z.wmW1 <- crop(precipw.z,east.z.wmW)
    elev.z.wmW1    <- crop(elev.z,east.z.wmW)
    north.z.wmW1   <- crop(north.z,east.z.wmW)
    rough.z.wmW1   <- crop(rough.z,east.z.wmW)
    wash.z.wmW1    <- crop(wash.z,east.z.wmW)
    bedrock.z.wmW1 <- crop(bedrock.z,east.z.wmW)
    vega.z.wmW1    <- crop(vega.z,east.z.wmW)
    roada.z.wmW1   <- crop(roada.z,east.z.wmW)
    precips.z.wmW <- mask(precips.z.wmW1,east.z.wmW)
    precipw.z.wmW <- mask(precipw.z.wmW1,east.z.wmW)
    elev.z.wmW <- mask(elev.z.wmW1,east.z.wmW)
    north.z.wmW <- mask(north.z.wmW1,east.z.wmW)
    rough.z.wmW <- mask(rough.z.wmW1,east.z.wmW)
    wash.z.wmW <- mask(wash.z.wmW1,east.z.wmW)
    bedrock.z.wmW <- mask(bedrock.z.wmW1,east.z.wmW)
    vega.z.wmW <- mask(vega.z.wmW1,east.z.wmW)
    roada.z.wmW <- mask(roada.z.wmW1,east.z.wmW)

    #Western Mojave - east half
    precips.z.wmE1 <- crop(precips.z,east.z.wmE)
    precipw.z.wmE1 <- crop(precipw.z,east.z.wmE)
    elev.z.wmE1    <- crop(elev.z,east.z.wmE)
    north.z.wmE1   <- crop(north.z,east.z.wmE)
    rough.z.wmE1   <- crop(rough.z,east.z.wmE)
    wash.z.wmE1    <- crop(wash.z,east.z.wmE)
    bedrock.z.wmE1 <- crop(bedrock.z,east.z.wmE)
    vega.z.wmE1    <- crop(vega.z,east.z.wmE)
    roada.z.wmE1   <- crop(roada.z,east.z.wmE)
    precips.z.wmE <- mask(precips.z.wmE1,east.z.wmE)
    precipw.z.wmE <- mask(precipw.z.wmE1,east.z.wmE)
    elev.z.wmE <- mask(elev.z.wmE1,east.z.wmE)
    north.z.wmE <- mask(north.z.wmE1,east.z.wmE)
    rough.z.wmE <- mask(rough.z.wmE1,east.z.wmE)
    wash.z.wmE <- mask(wash.z.wmE1,east.z.wmE)
    bedrock.z.wmE <- mask(bedrock.z.wmE1,east.z.wmE)
    vega.z.wmE <- mask(vega.z.wmE1,east.z.wmE)
    roada.z.wmE <- mask(roada.z.wmE1,east.z.wmE)
    
    #Eastern Mojave
    precips.z.em1 <- crop(precips.z,east.z.em)
    precipw.z.em1 <- crop(precipw.z,east.z.em)
    elev.z.em1    <- crop(elev.z,east.z.em)
    north.z.em1   <- crop(north.z,east.z.em)
    rough.z.em1   <- crop(rough.z,east.z.em)
    wash.z.em1    <- crop(wash.z,east.z.em)
    bedrock.z.em1 <- crop(bedrock.z,east.z.em)
    vega.z.em1    <- crop(vega.z,east.z.em)
    roada.z.em1   <- crop(roada.z,east.z.em)
    precips.z.em <- mask(precips.z.em1,east.z.em)
    precipw.z.em <- mask(precipw.z.em1,east.z.em)
    elev.z.em <- mask(elev.z.em1,east.z.em)
    north.z.em <- mask(north.z.em1,east.z.em)
    rough.z.em <- mask(rough.z.em1,east.z.em)
    wash.z.em <- mask(wash.z.em1,east.z.em)
    bedrock.z.em <- mask(bedrock.z.em1,east.z.em)
    vega.z.em <- mask(vega.z.em1,east.z.em)
    roada.z.em <- mask(roada.z.em1,east.z.em)
    
    #Northeastern Mojave
    precips.z.nm1 <- crop(precips.z,east.z.nm)
    precipw.z.nm1 <- crop(precipw.z,east.z.nm)
    elev.z.nm1    <- crop(elev.z,east.z.nm)
    north.z.nm1   <- crop(north.z,east.z.nm)
    rough.z.nm1   <- crop(rough.z,east.z.nm)
    wash.z.nm1    <- crop(wash.z,east.z.nm)
    bedrock.z.nm1 <- crop(bedrock.z,east.z.nm)
    vega.z.nm1    <- crop(vega.z,east.z.nm)
    roada.z.nm1   <- crop(roada.z,east.z.nm)
    precips.z.nm <- mask(precips.z.nm1,east.z.nm)
    precipw.z.nm <- mask(precipw.z.nm1,east.z.nm)
    elev.z.nm <- mask(elev.z.nm1,east.z.nm)
    north.z.nm <- mask(north.z.nm1,east.z.nm)
    rough.z.nm <- mask(rough.z.nm1,east.z.nm)
    wash.z.nm <- mask(wash.z.nm1,east.z.nm)
    bedrock.z.nm <- mask(bedrock.z.nm1,east.z.nm)
    vega.z.nm <- mask(vega.z.nm1,east.z.nm)
    roada.z.nm <- mask(roada.z.nm1,east.z.nm)

#-----------------------------------------------------------------------------------------------# 
# Create RasterStacks for predictions in each unit
#-----------------------------------------------------------------------------------------------#  
    
  #Putting covariates in the same order they appear in the model (cov_lam)
  cov.cdN <- stack(east.z.cdN,north.z.cdN,bedrock.z.cdN,elev.z.cdN,precips.z.cdN,
                   precipw.z.cdN,roada.z.cdN,rough.z.cdN,vega.z.cdN,wash.z.cdN)
  cov.cdS <- stack(east.z.cdS,north.z.cdS,bedrock.z.cdS,elev.z.cdS,precips.z.cdS,
                   precipw.z.cdS,roada.z.cdS,rough.z.cdS,vega.z.cdS,wash.z.cdS)
  cov.wmW <- stack(east.z.wmW,north.z.wmW,bedrock.z.wmW,elev.z.wmW,precips.z.wmW,
                   precipw.z.wmW,roada.z.wmW,rough.z.wmW,vega.z.wmW,wash.z.wmW)
  cov.wmE <- stack(east.z.wmE,north.z.wmE,bedrock.z.wmE,elev.z.wmE,precips.z.wmE,
                   precipw.z.wmE,roada.z.wmE,rough.z.wmE,vega.z.wmE,wash.z.wmE)
  cov.em <- stack(east.z.em,north.z.em,bedrock.z.em,elev.z.em,precips.z.em,
                   precipw.z.em,roada.z.em,rough.z.em,vega.z.em,wash.z.em)
  cov.nm <- stack(east.z.nm,north.z.nm,bedrock.z.nm,elev.z.nm,precips.z.nm,
                   precipw.z.nm,roada.z.nm,rough.z.nm,vega.z.nm,wash.z.nm)

  #Convert stacks to dataframes and then to matrices:
    
    #Note: each column is a layer, rows are in order of raster cells
    covm.cdN <- as.matrix(as.data.frame(cov.cdN)) 
    covm.cdS <- as.matrix(as.data.frame(cov.cdS))
    covm.wmW <- as.matrix(as.data.frame(cov.wmW))
    covm.wmE <- as.matrix(as.data.frame(cov.wmE))
    covm.em <- as.matrix(as.data.frame(cov.em))
    covm.nm <- as.matrix(as.data.frame(cov.nm))

  #Get indices of cell numbers where all covariates have values (no NAs) 
    
    cdN.ind <- which(apply(covm.cdN,1,function(x) sum(is.na(x))==0))   
    cdS.ind <- which(apply(covm.cdS,1,function(x) sum(is.na(x))==0))
    wmW.ind <- which(apply(covm.wmW,1,function(x) sum(is.na(x))==0))
    wmE.ind <- which(apply(covm.wmE,1,function(x) sum(is.na(x))==0))
    em.ind <- which(apply(covm.em,1,function(x) sum(is.na(x))==0))
    nm.ind <- which(apply(covm.nm,1,function(x) sum(is.na(x))==0))
    
  #Remove NA rows from data matrix
      
    covm.cdN1 <- covm.cdN[cdN.ind,]  #15261 cells
    covm.cdS1 <- covm.cdS[cdS.ind,]  #15462 
    covm.wmW1 <- covm.wmW[wmW.ind,]  #17148
    covm.wmE1 <- covm.wmE[wmE.ind,]  #19601
    covm.em1 <- covm.em[em.ind,]     #21373
    covm.nm1 <- covm.nm[nm.ind,]     #19986

#-----------------------------------------------------------------------------------------------# 
# Calculating predicted densities
#-----------------------------------------------------------------------------------------------#     
 
#Extracting posterior samples for regression coefficients in abundance model
  
  #Note: order of covariates in RasterStacks and in model parameters (cov_lam) must be the same
  #east, north, bedrock, elev, precips, precipw, road, rough, veg, wash
  
  iter <- seq(1,6000,by=6)
  
  #Regression coefficients in the abundance model (on log scale)
  betas <- as.matrix(posterior6000[iter,grep('beta_lam',colnames(posterior6000))])
  
  #Recovery unit intercepts (density, on real scale, in year 2000)
  mus <- as.matrix(posterior6000[iter,grep('mu_recov',colnames(posterior6000))])
  
  #Regression coefficient for yearly trends in the abundance model (on log scale)
  trends <- as.matrix(posterior6000[iter,grep('trend',colnames(posterior6000))])
  
  #Random effect of year[2018] (unique for each recovery unit)
  sd_recov <- as.matrix(posterior6000[iter,grep('sd_recov',colnames(posterior6000))])
  variates <- as.matrix(posterior6000[iter,grep('recov_l_raw',colnames(posterior6000))])
    #order: [1,1], [2,1], [3,1], [4,1], [1,2], ...[4,18]
  variates18 <- as.matrix(variates[,grep('18',colnames(variates))])
  re.cd.v <- as.vector(sd_recov*variates18[,1])
  re.cdN <- matrix(rep(re.cd.v,nrow(covm.cdN1)),byrow=TRUE,nrow=nrow(covm.cdN1))
  re.cdS <- matrix(rep(re.cd.v,nrow(covm.cdS1)),byrow=TRUE,nrow=nrow(covm.cdS1))
  re.wm.v <- as.vector(sd_recov*variates18[,4])
  re.wmW <- matrix(rep(re.wm.v,nrow(covm.wmW1)),byrow=TRUE,nrow=nrow(covm.wmW1))
  re.wmE <- matrix(rep(re.wm.v,nrow(covm.wmE1)),byrow=TRUE,nrow=nrow(covm.wmE1))
  re.em.v <- as.vector(sd_recov*variates18[,2])
  re.em <- matrix(rep(re.em.v,nrow(covm.em1)),byrow=TRUE,nrow=nrow(covm.em1))
  re.nm.v <- as.vector(sd_recov*variates18[,3])
  re.nm <- matrix(rep(re.nm.v,nrow(covm.nm1)),byrow=TRUE,nrow=nrow(covm.nm1))
  
#Calculating predicted densities
  #log(pred.yr18) = log(mus[RU]) + 18 %*% trend[RU] + covm %*% betas + randomeffect[RU,18] 
  #order of RUs: CD, EM, NM, WM
  
  #CD - North
    cd.mu <- as.matrix(log(mus[,1]))
    cd.trend <- as.matrix(trends[,1])
    cdN.predl <- matrix(1,nrow=nrow(covm.cdN1),ncol=1) %*% t(cd.mu) + matrix(18,nrow=nrow(covm.cdN1),ncol=1) %*% t(cd.trend) +
                 covm.cdN1 %*% t(betas) + re.cdN
    cdN.pred <- exp(cdN.predl)
    cdN.predm <- apply(cdN.pred,1,median)
    cdN.preds <- east.z.cdN
    cdN.preds[] <- NA
    names(cdN.preds) <- 'preds'
    cdN.preds[cdN.ind] <- cdN.predm
    #plot(cdN.preds)
  
  #CD - South
    cdS.predl <- matrix(1,nrow=nrow(covm.cdS1),ncol=1) %*% t(cd.mu) + matrix(18,nrow=nrow(covm.cdS1),ncol=1) %*% t(cd.trend) +
                 covm.cdS1 %*% t(betas) + re.cdS
    cdS.pred <- exp(cdS.predl)
    cdS.predm <- apply(cdS.pred,1,median)
    cdS.preds <- east.z.cdS
    cdS.preds[] <- NA
    names(cdS.preds) <- 'preds'
    cdS.preds[cdS.ind] <- cdS.predm
    #plot(cdS.preds)
    
  #Western Mojave - West    
    wm.mu <- as.matrix(log(mus[,4]))
    wm.trend <- as.matrix(trends[,4])
    wmW.predl <- matrix(1,nrow=nrow(covm.wmW1),ncol=1) %*% t(wm.mu) + matrix(18,nrow=nrow(covm.wmW1),ncol=1) %*% t(wm.trend) +
                 covm.wmW1 %*% t(betas) + re.wmW
    wmW.pred <- exp(wmW.predl)
    wmW.predm <- apply(wmW.pred,1,median)
    wmW.preds <- east.z.wmW
    wmW.preds[] <- NA
    names(wmW.preds) <- 'preds'
    wmW.preds[wmW.ind] <- wmW.predm
    #plot(wmW.preds)
    
  #Western Mojave - East    
    wmE.predl <- matrix(1,nrow=nrow(covm.wmE1),ncol=1) %*% t(wm.mu) + matrix(18,nrow=nrow(covm.wmE1),ncol=1) %*% t(wm.trend) +
                 covm.wmE1 %*% t(betas) + re.wmE
    wmE.pred <- exp(wmE.predl)
    wmE.predm <- apply(wmE.pred,1,median)
    wmE.preds <- east.z.wmE
    wmE.preds[] <- NA
    names(wmE.preds) <- 'preds'
    wmE.preds[wmE.ind] <- wmE.predm
    #plot(wmE.preds)    
    
  #Eastern Mojave
    em.mu <- as.matrix(log(mus[,2]))
    em.trend <- as.matrix(trends[,2])
    em.predl <- matrix(1,nrow=nrow(covm.em1),ncol=1) %*% t(em.mu) + matrix(18,nrow=nrow(covm.em1),ncol=1) %*% t(em.trend) +
                covm.em1 %*% t(betas) + re.em
    em.pred <- exp(em.predl)
    em.predm <- apply(em.pred,1,median)
    em.preds <- east.z.em
    em.preds[] <- NA
    names(em.preds) <- 'preds'
    em.preds[em.ind] <- em.predm
    #plot(em.preds)    
    
  #Northeastern Mojave
    nm.mu <- as.matrix(log(mus[,3]))
    nm.trend <- as.matrix(trends[,3])
    nm.predl <- matrix(1,nrow=nrow(covm.nm1),ncol=1) %*% t(nm.mu) + matrix(18,nrow=nrow(covm.nm1),ncol=1) %*% t(nm.trend) +
                covm.nm1 %*% t(betas) + re.nm
    nm.pred <- exp(nm.predl)
    nm.predm <- apply(nm.pred,1,median)
    nm.preds <- east.z.nm
    nm.preds[] <- NA
    names(nm.preds) <- 'preds'
    nm.preds[nm.ind] <- nm.predm
    #plot(nm.preds) 

 #Merge rasters and plot   
    r6 <- list(cdN.preds,cdS.preds,wmW.preds,wmE.preds,em.preds,nm.preds)  
    pred6 <- do.call(merge,r6)
    
    starts <- surveys[,c('segmentID','start_e','start_n')]
    startlocs <- SpatialPointsDataFrame(coords=starts[,2:3],proj4string=crs(east),data=starts)

    #jpeg('Figures/PredictedDensities_wREs.jpg',width=6.5,height=6.5,units='in',res=600)
      par(mar=c(1,1,1,1))
      plot(pred6)
      # plot(roads.nv,col='gray50',add=TRUE)
      # plot(roads.ca,col='gray50',add=TRUE)
      # plot(roads.az,col='gray50',add=TRUE)
      # plot(roads.ut,col='gray50',add=TRUE)
      plot(states,lwd=2,add=TRUE)
      plot(rus,lwd=3,add=TRUE)
      # points(startlocs,pch=19,cex=0.3,col='dodgerblue')
    #dev.off() 
