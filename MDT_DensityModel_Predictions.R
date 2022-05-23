#===============================================================================================# 

 ## Authors: ER Zylstra, RJ Steidl

 ## Project: MDT line-distance sampling - estimating density (spatial variation and trends)

  # Model structure adapted from Chelgren et al. 2011, Herpetol. Conserv. Biol. 6:175-190
  # and Shirk et al. 2014, Diversity and Distributions 20:1186-1199
  
  ## Relevant files: 
  
    # 1. MDT_Import_Format_Data.R -- import all the files needed for analysis, format data, and create MDT_Data.Rdata
    # 2. MDT_DensityModel.R -- run model in STAN, summarize results, produce figures (study area map, trends)
    # 3. MDT_DensityModel.stan -- STAN model
    # 4. MDT_DensityModel_Predictions.R -- predict density of tortoises across their range (produce heat map)

#===============================================================================================#  

#-----------------------------------------------------------------------------------------------# 
# Set working directory, load packages
#-----------------------------------------------------------------------------------------------# 

  # Set working directory

    # setwd()

  # Load packages

    library(raster)
    library(rgdal)
    library(ggplot2)
    library(dplyr)

  # rm(list=ls())
  
#-----------------------------------------------------------------------------------------------# 
# Import data
#-----------------------------------------------------------------------------------------------# 
  
  # Import survey and associated covariate data (so I can standardize the entire raster appropriately)

    load('MDT_Data.Rdata')
    
  # Output from Stan model (and associated covariate data that went into the model)
    
    load('MDT_Density_2001-2020.Rdata')  
    posterior6000 <- as.matrix(out)
    
  # Covariate data
    
    precip.s <- raster('Covariates/30yr_may_oct_prcp_mm.tif')            #Summer precip (May-Oct, mm), 30yr norms
    precip.w <- raster('Covariates/30yr_nov_apr_prcp_mm.tif')            #Winter precip (Nov-Apr, mm), 30yr norms
    elev <- raster('Covariates/avg_elevation.tif')                       #Elevation (m)
    slope <- raster('Covariates/avg_slope.tif')                          #Slope (degrees)
    north <- raster('Covariates/avg_northness.tif')                      #Aspect:northness 
    east <- raster('Covariates/avg_eastness.tif')                        #Aspect:eastness 
    rough <- raster('Covariates/average_surface_roughness_snapped.tif')  #Average surface roughness
    wash <- raster('Covariates/wash_proportion.tif')                     #Wash density (proportion of each pixel predicted to be "wash")
    bedrock <- raster('Covariates/bedrock_depth.tif')                    #Depth to bedrock (cm)
    veg.p <- raster('Covariates/perennial_vegetation.tif')               #Perennial plant cover: NDVI values in 2019 (15 Jun - 31 Aug), very dry year 
    veg.a <- raster('Covariates/annProx_snapped.tif')                    #Annual plant growth potential (MODIS, compare wet(2005)-dry(2002) years)
    road.a <- raster('Covariates/dist_road.tif')                         #Distance to any road (m)
    
  # Shapefile: Recovery units
    
    rus <- readOGR(dsn='Covariates/Revised Recovery Units',layer='2011RecoveryUnits')
    rus <- spTransform(rus,crs(east))
    
  # Shapefile: TCAs
    
    tcas <- readOGR(dsn='Covariates/TCAs',layer='All_Strata')
    tcas <- spTransform(tcas,crs(east))
    
  # Shapefile: State boundaries
    
    states <- readOGR(dsn='Covariates/US_states_GIS',layer='cb_2017_us_state_500k')
    states <- spTransform(states,crs(east))

  # Shapefile: Prediction Area 
    # Includes any cells within 4 RUs that meet covariate criteria below
    # Excludes 1 sq km cells where one or more covariates values are NA
    # Excludes cells where any covariate is more than 10% outside the range of 
      # values observed at survey locations
    
    predarea <- readOGR(dsn='Covariates/PredictionArea_Full',layer='PredictionArea_Full')
    predarea <- spTransform(predarea,crs(east))    

  # Impervious surfaces
    
    # Proportion impervious for each 1 sq km cell in 2001
    imp01 <- raster('Covariates/ImpervSurfaces/impervious_raw_2001.tif')   
    
    # Proportion impervious for each 1 sq km cell in 2019
    imp19 <- raster('Covariates/ImpervSurfaces/impervious_raw_2019.tif')   

#-----------------------------------------------------------------------------------------------# 
# Log covariate rasters (where needed) and impute missing values
#-----------------------------------------------------------------------------------------------# 

  # Wash layer is missing some values in the north-central part of the study area
  # Will impute the mean value, like we did before extracting values for survey locations
    
    wash.missing <- wash
    wash[is.na(wash.missing)] <- mean(values(wash.missing),na.rm=TRUE)  
    # Adjusting minimum value of wash, which was very small, but >0 for all surveys
    wash[wash < min(surveys$wash)] <- min(surveys$wash)
      
  # Log covariates that were logged in the density model
    
    rough.log <- log(rough)
    wash.log <- log(wash)
    bedrock.log <- log(bedrock)
    road.a.log <- log(road.a)    

#-----------------------------------------------------------------------------------------------# 
# Standardize values in covariate rasters
#-----------------------------------------------------------------------------------------------# 

  # Covariates that were log-transformed before standardizing: rough, wash, bedrock, road.a
    
    surveys$rough.log <- log(surveys$rough)
    surveys$wash.log <- log(surveys$wash)
    surveys$bedrock.log <- log(surveys$bedrock)
    surveys$road.a.log <- log(surveys$road.a)
    
  # Getting means, SDs for each covariate in the density model
    
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

  # Standardize rasters 
    
    precips.z <- (precip.s - precips.mn)/precips.sd
    precipw.z <- (precip.w - precipw.mn)/precipw.sd
    elev.z <- (elev - elev.mn)/elev.sd
    rough.z <- (rough.log - rough.l.mn)/rough.l.sd
    east.z <- (east - east.mn)/east.sd
    north.z <- (north - north.mn)/north.sd
    wash.z <- (wash.log - wash.l.mn)/wash.l.sd
    bedrock.z <- (bedrock.log - bedrock.l.mn)/bedrock.l.sd
    vega.z <- (veg.a - vega.mn)/vega.sd
    roada.z <- (road.a.log - roada.l.mn)/roada.l.sd      

#-----------------------------------------------------------------------------------------------#     
# Crop standardized covariate rasters to prediction areas within each recovery unit 
#-----------------------------------------------------------------------------------------------# 

  #Polygon for each recovery unit
    
    cd <- subset(rus,Unit_Name=='Colorado Desert')
    em <- subset(rus,Unit_Name=='Eastern Mojave')
    nm <- subset(rus,Unit_Name=='Northeastern Mojave')
    wm <- subset(rus,Unit_Name=='Western Mojave')
      
  #Create shapefile with prediction area in each recovery unit
    
    pred.cd <- crop(predarea,cd)
    pred.em <- crop(predarea,em)
    pred.nm <- crop(predarea,nm)
    pred.wm <- crop(predarea,wm)

  #Crop (and mask) all covariate layers in CD
    
    precips.cd <- crop(precips.z,pred.cd)
    precips.cd <- mask(precips.cd,pred.cd)
    precipw.cd <- crop(precipw.z,pred.cd)
    precipw.cd <- mask(precipw.cd,pred.cd)
    elev.cd <- crop(elev.z,pred.cd)
    elev.cd <- mask(elev.cd,pred.cd)
    north.cd <- crop(north.z,pred.cd)
    north.cd <- mask(north.cd,pred.cd)
    east.cd <- crop(east.z,pred.cd)
    east.cd <- mask(east.cd,pred.cd)
    rough.cd <- crop(rough.z,pred.cd)
    rough.cd <- mask(rough.cd,pred.cd)
    wash.cd <- crop(wash.z,pred.cd)
    wash.cd <- mask(wash.cd,pred.cd)
    bedrock.cd <- crop(bedrock.z,pred.cd)
    bedrock.cd <- mask(bedrock.cd,pred.cd)
    vega.cd <- crop(vega.z,pred.cd)
    vega.cd <- mask(vega.cd,pred.cd)
    roada.cd <- crop(roada.z,pred.cd)
    roada.cd <- mask(roada.cd,pred.cd)
    
    # Crop (and mask) impervious surface layers in CD
    imp01.cd <- crop(imp01,pred.cd)
    imp01.cd <- mask(imp01.cd,pred.cd)
    imp19.cd <- crop(imp19,pred.cd)
    imp19.cd <- mask(imp19.cd,pred.cd)

  # Crop (and mask) all covariate layers in EM
    
    precips.em <- crop(precips.z,pred.em)
    precips.em <- mask(precips.em,pred.em)
    precipw.em <- crop(precipw.z,pred.em)
    precipw.em <- mask(precipw.em,pred.em)
    elev.em <- crop(elev.z,pred.em)
    elev.em <- mask(elev.em,pred.em)
    north.em <- crop(north.z,pred.em)
    north.em <- mask(north.em,pred.em)
    east.em <- crop(east.z,pred.em)
    east.em <- mask(east.em,pred.em)
    rough.em <- crop(rough.z,pred.em)
    rough.em <- mask(rough.em,pred.em)
    wash.em <- crop(wash.z,pred.em)
    wash.em <- mask(wash.em,pred.em)
    bedrock.em <- crop(bedrock.z,pred.em)
    bedrock.em <- mask(bedrock.em,pred.em)
    vega.em <- crop(vega.z,pred.em)
    vega.em <- mask(vega.em,pred.em)
    roada.em <- crop(roada.z,pred.em)
    roada.em <- mask(roada.em,pred.em)
    
    # Crop (and mask) impervious surface layers in EM
    imp01.em <- crop(imp01,pred.em)
    imp01.em <- mask(imp01.em,pred.em)
    imp19.em <- crop(imp19,pred.em)
    imp19.em <- mask(imp19.em,pred.em)
    
  # Crop (and mask) all covariate layers in NM
    
    precips.nm <- crop(precips.z,pred.nm)
    precips.nm <- mask(precips.nm,pred.nm)
    precipw.nm <- crop(precipw.z,pred.nm)
    precipw.nm <- mask(precipw.nm,pred.nm)
    elev.nm <- crop(elev.z,pred.nm)
    elev.nm <- mask(elev.nm,pred.nm)
    north.nm <- crop(north.z,pred.nm)
    north.nm <- mask(north.nm,pred.nm)
    east.nm <- crop(east.z,pred.nm)
    east.nm <- mask(east.nm,pred.nm)
    rough.nm <- crop(rough.z,pred.nm)
    rough.nm <- mask(rough.nm,pred.nm)
    wash.nm <- crop(wash.z,pred.nm)
    wash.nm <- mask(wash.nm,pred.nm)
    bedrock.nm <- crop(bedrock.z,pred.nm)
    bedrock.nm <- mask(bedrock.nm,pred.nm)
    vega.nm <- crop(vega.z,pred.nm)
    vega.nm <- mask(vega.nm,pred.nm)
    roada.nm <- crop(roada.z,pred.nm)
    roada.nm <- mask(roada.nm,pred.nm)
    
    # Crop (and mask) impervious surface layers in NM
    imp01.nm <- crop(imp01,pred.nm)
    imp01.nm <- mask(imp01.nm,pred.nm)
    imp19.nm <- crop(imp19,pred.nm)
    imp19.nm <- mask(imp19.nm,pred.nm)
    
  # Crop (and mask) all covariate layers in WM
    
    precips.wm <- crop(precips.z,pred.wm)
    precips.wm <- mask(precips.wm,pred.wm)
    precipw.wm <- crop(precipw.z,pred.wm)
    precipw.wm <- mask(precipw.wm,pred.wm)
    elev.wm <- crop(elev.z,pred.wm)
    elev.wm <- mask(elev.wm,pred.wm)
    north.wm <- crop(north.z,pred.wm)
    north.wm <- mask(north.wm,pred.wm)
    east.wm <- crop(east.z,pred.wm)
    east.wm <- mask(east.wm,pred.wm)
    rough.wm <- crop(rough.z,pred.wm)
    rough.wm <- mask(rough.wm,pred.wm)
    wash.wm <- crop(wash.z,pred.wm)
    wash.wm <- mask(wash.wm,pred.wm)
    bedrock.wm <- crop(bedrock.z,pred.wm)
    bedrock.wm <- mask(bedrock.wm,pred.wm)
    vega.wm <- crop(vega.z,pred.wm)
    vega.wm <- mask(vega.wm,pred.wm)
    roada.wm <- crop(roada.z,pred.wm)
    roada.wm <- mask(roada.wm,pred.wm) 
  
    # Crop (and mask) impervious surface layers in WM
    imp01.wm <- crop(imp01,pred.wm)
    imp01.wm <- mask(imp01.wm,pred.wm)
    imp19.wm <- crop(imp19,pred.wm)
    imp19.wm <- mask(imp19.wm,pred.wm)

#-----------------------------------------------------------------------------------------------# 
# Define (binary) impervious surfaces
#-----------------------------------------------------------------------------------------------#   

  is.threshold <- 0.4

  # Get indices of cell numbers where impervious surface values are below threshold
  
    imp01.cd.ind <- which(values(imp01.cd) < is.threshold)
    imp19.cd.ind <- which(values(imp19.cd) < is.threshold)
    imp01.em.ind <- which(values(imp01.em) < is.threshold)
    imp19.em.ind <- which(values(imp19.em) < is.threshold)
    imp01.nm.ind <- which(values(imp01.nm) < is.threshold)
    imp19.nm.ind <- which(values(imp19.nm) < is.threshold)
    imp01.wm.ind <- which(values(imp01.wm) < is.threshold)
    imp19.wm.ind <- which(values(imp19.wm) < is.threshold)

  # Because impervious surface layer is already cropped/masked by prediction area, 
  # we can use the indexes of cells where impervious surface value < threshold to make predictions    
        
#-----------------------------------------------------------------------------------------------# 
# Organize covariates for prediction in each grid cell
#-----------------------------------------------------------------------------------------------#  
    
  # Putting covariates in the same order they appear in the model (cov_lam)
    
    covs.cd <- stack(east.cd,north.cd,bedrock.cd,elev.cd,precips.cd,precipw.cd,roada.cd,rough.cd,vega.cd,wash.cd)
    covs.em <- stack(east.em,north.em,bedrock.em,elev.em,precips.em,precipw.em,roada.em,rough.em,vega.em,wash.em)
    covs.nm <- stack(east.nm,north.nm,bedrock.nm,elev.nm,precips.nm,precipw.nm,roada.nm,rough.nm,vega.nm,wash.nm)
    covs.wm <- stack(east.wm,north.wm,bedrock.wm,elev.wm,precips.wm,precipw.wm,roada.wm,rough.wm,vega.wm,wash.wm)
    
  # Convert each stack to a dataframe and then to a matrix:
    # Note: each column is a layer, rows are in order of raster cells
    
    covm.cd <- as.matrix(as.data.frame(covs.cd)) 
    covm.em <- as.matrix(as.data.frame(covs.em)) 
    covm.nm <- as.matrix(as.data.frame(covs.nm)) 
    covm.wm <- as.matrix(as.data.frame(covs.wm)) 

  # Check that for a given raster cell, no covariate values are NA or all are NA
  # (b/c that's how we defined the prediction area used to crop covariate rasters
    
    table(apply(covm.cd,1,function(x) sum(is.na(x))))
    table(apply(covm.em,1,function(x) sum(is.na(x))))
    table(apply(covm.nm,1,function(x) sum(is.na(x))))
    table(apply(covm.wm,1,function(x) sum(is.na(x))))
 
  #Remove rows from data matrix where covariates = NA or it's impervious
    
    covm01.cd <- covm.cd[imp01.cd.ind,]
    covm19.cd <- covm.cd[imp19.cd.ind,]
    covm01.em <- covm.em[imp01.em.ind,]
    covm19.em <- covm.em[imp19.em.ind,]
    covm01.nm <- covm.nm[imp01.nm.ind,]
    covm19.nm <- covm.nm[imp19.nm.ind,]
    covm01.wm <- covm.wm[imp01.wm.ind,]
    covm19.wm <- covm.wm[imp19.wm.ind,]

#-----------------------------------------------------------------------------------------------# 
# Calculating predicted densities for each raster cell in 2001 (after removing impervious surfaces)
#-----------------------------------------------------------------------------------------------#     

  # Note: Not including random yearly effects around range-wide trend, which are more likely to 
    # capture differences in areas sampled each year than true density changes

  # Extracting posterior samples for regression coefficients in abundance model
    # Note: order of covariates in RasterStacks and in model parameters (cov_lam) must be the same
    # east, north, bedrock, elev, precips, precipw, road, rough, veg, wash

    iter <- seq(1,6000,by=6)
    
    # Regression coefficients in the abundance model (on log scale)
    
      betas <- as.matrix(posterior6000[iter,grep('beta_lam',colnames(posterior6000))])
    
    #Recovery unit intercepts (density on real scale in year 2000)
    
      mus <- as.matrix(posterior6000[iter,grep('mu_recov',colnames(posterior6000))])
    
    #Regression coefficient for yearly trends in the abundance model (on log scale)
      
      trends <- as.matrix(posterior6000[iter,grep('trend',colnames(posterior6000))])

  # Calculating predicted densities in 2001
    # log(pred.yr1) = log(mus[RU]) + 1 %*% trend[RU] + covm01.ru %*% betas
    # order of RUs: CD, EM, NM, WM
    
    # CD
      
      cd.mu <- as.matrix(log(mus[,1]))
      cd.trend <- as.matrix(trends[,1])
      cd.predl01 <- matrix(1,nrow=nrow(covm01.cd),ncol=1) %*% t(cd.mu) + 
                    matrix(1,nrow=nrow(covm01.cd),ncol=1) %*% t(cd.trend) +
                    covm01.cd %*% t(betas)
      cd.pred01 <- exp(cd.predl01)
      #Median expected density in each cell
      cd.predm01 <- apply(cd.pred01,1,median) 
      #Creating raster with median density values
      cd.preds01 <- east.cd
      cd.preds01[] <- NA
      cd.preds01[imp01.cd.ind] <- cd.predm01

    # EM
      
      em.mu <- as.matrix(log(mus[,2]))
      em.trend <- as.matrix(trends[,2]) 
      em.predl01 <- matrix(1,nrow=nrow(covm01.em),ncol=1) %*% t(em.mu) + 
                    matrix(1,nrow=nrow(covm01.em),ncol=1) %*% t(em.trend) +
                    covm01.em %*% t(betas)
      em.pred01 <- exp(em.predl01)
      em.predm01 <- apply(em.pred01,1,median)
      em.preds01 <- east.em
      em.preds01[] <- NA
      em.preds01[imp01.em.ind] <- em.predm01

    # NM
      
      nm.mu <- as.matrix(log(mus[,3]))
      nm.trend <- as.matrix(trends[,3])      
      nm.predl01 <- matrix(1,nrow=nrow(covm01.nm),ncol=1) %*% t(nm.mu) + 
                    matrix(1,nrow=nrow(covm01.nm),ncol=1) %*% t(nm.trend) +
                    covm01.nm %*% t(betas)
      nm.pred01 <- exp(nm.predl01)
      nm.predm01 <- apply(nm.pred01,1,median)
      nm.preds01 <- east.nm
      nm.preds01[] <- NA
      nm.preds01[imp01.nm.ind] <- nm.predm01

    # WM
      
      wm.mu <- as.matrix(log(mus[,4]))
      wm.trend <- as.matrix(trends[,4])      
      wm.predl01 <- matrix(1,nrow=nrow(covm01.wm),ncol=1) %*% t(wm.mu) + 
                    matrix(1,nrow=nrow(covm01.wm),ncol=1) %*% t(wm.trend) +
                    covm01.wm %*% t(betas)
      wm.pred01 <- exp(wm.predl01)
      wm.predm01 <- apply(wm.pred01,1,median)
      wm.preds01 <- east.wm
      wm.preds01[] <- NA
      wm.preds01[imp01.wm.ind] <- wm.predm01

  # Raster for all recovery units    
      
    r01 <- list(cd.preds01, em.preds01, nm.preds01, wm.preds01)  
    preds01 <- do.call(merge,r01) 
    plot(preds01)

#-----------------------------------------------------------------------------------------------# 
# Calculating predicted densities for each raster cell in 2020 (after removing 2019 impervious surfaces)
#-----------------------------------------------------------------------------------------------#     
  
  # Calculating predicted densities in 2020 
    # log(pred.yr20) = log(mus[RU]) + 20 %*% trend[RU] + covm19.ru %*% betas
    # order of RUs: CD, EM, NM, WM
    
    # CD
    
      cd.predl20 <- matrix(1,nrow=nrow(covm19.cd),ncol=1) %*% t(cd.mu) + 
                    matrix(20,nrow=nrow(covm19.cd),ncol=1) %*% t(cd.trend) +
                    covm19.cd %*% t(betas)
      cd.pred20 <- exp(cd.predl20)
      cd.predm20 <- apply(cd.pred20,1,median)
      cd.preds20 <- east.cd
      cd.preds20[] <- NA
      cd.preds20[imp19.cd.ind] <- cd.predm20

    # EM
      
      em.predl20 <- matrix(1,nrow=nrow(covm19.em),ncol=1) %*% t(em.mu) + 
                    matrix(20,nrow=nrow(covm19.em),ncol=1) %*% t(em.trend) +
                    covm19.em %*% t(betas)
      em.pred20 <- exp(em.predl20)
      em.predm20 <- apply(em.pred20,1,median)
      em.preds20 <- east.em
      em.preds20[] <- NA
      em.preds20[imp19.em.ind] <- em.predm20

    # NM
      
      nm.predl20 <- matrix(1,nrow=nrow(covm19.nm),ncol=1) %*% t(nm.mu) + 
                    matrix(20,nrow=nrow(covm19.nm),ncol=1) %*% t(nm.trend) +
                    covm19.nm %*% t(betas)
      nm.pred20 <- exp(nm.predl20)
      nm.predm20 <- apply(nm.pred20,1,median)
      nm.preds20 <- east.nm
      nm.preds20[] <- NA
      nm.preds20[imp19.nm.ind] <- nm.predm20

    #WM
      
      wm.predl20 <- matrix(1,nrow=nrow(covm19.wm),ncol=1) %*% t(wm.mu) + 
                    matrix(20,nrow=nrow(covm19.wm),ncol=1) %*% t(wm.trend) +
                    covm19.wm %*% t(betas)
      wm.pred20 <- exp(wm.predl20)
      wm.predm20 <- apply(wm.pred20,1,median)
      wm.preds20 <- east.wm
      wm.preds20[] <- NA
      wm.preds20[imp19.wm.ind] <- wm.predm20

  # Raster for all recovery units 
      
    r20 <- list(cd.preds20, em.preds20, nm.preds20, wm.preds20)
    preds20 <- do.call(merge,r20)      
    plot(preds20)
    
  # Mean predicted density across all units
    
    all.pred20 <- rbind(cd.pred20, em.pred20, nm.pred20, wm.pred20)
    all.predm20 <- apply(all.pred20, 1, median)
    mean(all.predm20); quantile(all.predm20, probs = c(0.025, 0.975))
  
  # Calculating range-wide trend
  # Weighted average of trend estimates (weighted by the number of cells in each RU in 2020)
    
    n.cd <- nrow(cd.pred20)
    n.em <- nrow(em.pred20)
    n.nm <- nrow(nm.pred20)
    n.wm <- nrow(wm.pred20)
    n.tot <- n.cd + n.em + n.nm + n.wm
    otrend1000 <- apply(trends, 1, function(x) (n.cd*x[1] + n.em*x[2] + n.nm*x[3] + n.wm*x[4])/n.tot) 
    mean(otrend1000); quantile(otrend1000, probs = c(0.025, 0.975))
    
#-----------------------------------------------------------------------------------------------# 
# Summarizing 2001 and 2020 abundance and density in each recovery unit
#-----------------------------------------------------------------------------------------------#                

# Note: Estimated areas and abundances are for cells within the prediction area
# Excludes impervious surfaces (assuming 2019 values apply to 2020)
# Excludes cells with covariates values that are missing or significantly different than surveyed areas 
  
  # CD - 2001
    
    # Area
    area.cd01 <- length(imp01.cd.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.cd01 <- apply(cd.pred01,2,sum)
    A.cd01 <- round(median(A1000.cd01))
    A.cd01.sd <- round(sd(A1000.cd01))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.cd01 <- apply(cd.pred01,2,mean)
    D.cd01 <- round(median(D1000.cd01),2)
    D.cd01.sd <- round(sd(D1000.cd01),2)

  #CD - 2020
    
    # Area
    area.cd20 <- length(imp19.cd.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.cd20 <- apply(cd.pred20,2,sum)
    A.cd20 <- round(median(A1000.cd20))
    A.cd20.sd <- round(sd(A1000.cd20))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.cd20 <- apply(cd.pred20,2,mean)
    D.cd20 <- round(median(D1000.cd20),2)
    D.cd20.sd <- round(sd(D1000.cd20),2)
    
  # CD - Difference (2020-2001)
    
    diff1000.cd <- A1000.cd20 - A1000.cd01
    diff.cd <- round(median(diff1000.cd))
    diff.cd.sd <- round(sd(diff1000.cd))

  # EM - 2001
    
    # Area
    area.em01 <- length(imp01.em.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.em01 <- apply(em.pred01,2,sum)
    A.em01 <- round(median(A1000.em01))
    A.em01.sd <- round(sd(A1000.em01))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.em01 <- apply(em.pred01,2,mean)
    D.em01 <- round(median(D1000.em01),2)
    D.em01.sd <- round(sd(D1000.em01),2)

  #EM - 2020
    
    # Area
    area.em20 <- length(imp19.em.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.em20 <- apply(em.pred20,2,sum)
    A.em20 <- round(median(A1000.em20))
    A.em20.sd <- round(sd(A1000.em20))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.em20 <- apply(em.pred20,2,mean)
    D.em20 <- round(median(D1000.em20),2)
    D.em20.sd <- round(sd(D1000.em20),2)    
    
  # EM - Difference (2020-2001)
    
    diff1000.em <- A1000.em20 - A1000.em01
    diff.em <- round(median(diff1000.em))
    diff.em.sd <- round(sd(diff1000.em))
    
  # NM - 2001
    
    # Area
    area.nm01 <- length(imp01.nm.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.nm01 <- apply(nm.pred01,2,sum)
    A.nm01 <- round(median(A1000.nm01))
    A.nm01.sd <- round(sd(A1000.nm01))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.nm01 <- apply(nm.pred01,2,mean)
    D.nm01 <- round(median(D1000.nm01),2)
    D.nm01.sd <- round(sd(D1000.nm01),2)

  # NM - 2020
    
    # Area
    area.nm20 <- length(imp19.nm.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.nm20 <- apply(nm.pred20,2,sum)
    A.nm20 <- round(median(A1000.nm20))
    A.nm20.sd <- round(sd(A1000.nm20))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.nm20 <- apply(nm.pred20,2,mean)
    D.nm20 <- round(median(D1000.nm20),2)
    D.nm20.sd <- round(sd(D1000.nm20),2)    
    
  # NM - Difference (2020-2001)
    
    diff1000.nm <- A1000.nm20 - A1000.nm01
    diff.nm <- round(median(diff1000.nm))
    diff.nm.sd <- round(sd(diff1000.nm))
    
  # WM - 2001
    
    # Area
    area.wm01 <- length(imp01.wm.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.wm01 <- apply(wm.pred01,2,sum)
    A.wm01 <- round(median(A1000.wm01))
    A.wm01.sd <- round(sd(A1000.wm01))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.wm01 <- apply(wm.pred01,2,mean)
    D.wm01 <- round(median(D1000.wm01),2)
    D.wm01.sd <- round(sd(D1000.wm01),2)

  # WM - 2020
    
    # Area
    area.wm20 <- length(imp19.wm.ind)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.wm20 <- apply(wm.pred20,2,sum)
    A.wm20 <- round(median(A1000.wm20))
    A.wm20.sd <- round(sd(A1000.wm20))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.wm20 <- apply(wm.pred20,2,mean)
    D.wm20 <- round(median(D1000.wm20),2)
    D.wm20.sd <- round(sd(D1000.wm20),2)   
    
  # WM - Difference (2020-2001)
    
    diff1000.wm <- A1000.wm20 - A1000.wm01
    diff.wm <- round(median(diff1000.wm))
    diff.wm.sd <- round(sd(diff1000.wm))    

  # Across RUs - 2001
    
    all.pred01 <- rbind(cd.pred01, em.pred01, nm.pred01, wm.pred01)
    area.01 <- nrow(all.pred01)
    A1000.01 <- apply(all.pred01,2,sum)
    A.01 <- round(median(A1000.01))
    A.01.sd <- round(sd(A1000.01)) 
    D1000.01 <- apply(all.pred01,2,mean)
    D.01 <- round(median(D1000.01),2)
    D.01.sd <- round(sd(D1000.01),2)
    
  # Across RUs - 2020
    
    all.pred20 <- rbind(cd.pred20, em.pred20, nm.pred20, wm.pred20)
    area.20 <- nrow(all.pred20)
    A1000.20 <- apply(all.pred20,2,sum)
    A.20 <- round(median(A1000.20))
    A.20.sd <- round(sd(A1000.20)) 
    D1000.20 <- apply(all.pred20,2,mean)
    D.20 <- round(median(D1000.20),2)
    D.20.sd <- round(sd(D1000.20),2)    
    
  # Across RUs - Difference (2020-2001)
    
    diff1000 <- A1000.20 - A1000.01
    diff <- round(median(diff1000))
    diff.sd <- round(sd(diff1000))     

  abund <- data.frame(ru = c('CD', 'EM', 'NM', 'WM', 'Total'),
                      area.01 = c(area.cd01, area.em01, area.nm01, area.wm01, area.01), 
                      N.01 = c(A.cd01, A.em01, A.nm01, A.wm01, A.01),
                      N.01.sd = c(A.cd01.sd, A.em01.sd, A.nm01.sd, A.wm01.sd, A.01.sd),
                      area.20 = c(area.cd20, area.em20, area.nm20, area.wm20, area.20), 
                      N.20 = c(A.cd20, A.em20, A.nm20, A.wm20, A.20),
                      N.20.sd = c(A.cd20.sd, A.em20.sd, A.nm20.sd, A.wm20.sd, A.20.sd),
                      diff = c(diff.cd, diff.em, diff.nm, diff.wm, diff),
                      diff.sd = c(diff.cd.sd, diff.em.sd, diff.nm.sd, diff.wm.sd, diff.sd))
  abund

#-----------------------------------------------------------------------------------------------# 
# Calculating density and abundance in each TCA (prediction areas only)
#-----------------------------------------------------------------------------------------------#  

  # For simplicity, calculating the mean over median expected densities in each cell
  
    strata <- tcas@data$stratum[tcas@data$stratum!='RC']
    ag <- subset(tcas,stratum=='AG')
    jt <- subset(tcas,stratum=='JT')
    bd <- subset(tcas,stratum=='BD')
    ck <- subset(tcas,stratum=='CK')
    cm <- subset(tcas,stratum=='CM')
    cs <- subset(tcas,stratum=='CS')
    fe <- subset(tcas,stratum=='FE')
    fk <- subset(tcas,stratum=='FK')
    gb <- subset(tcas,stratum=='GB')
    iv <- subset(tcas,stratum=='IV')
    mm <- subset(tcas,stratum=='MM')
    or <- subset(tcas,stratum=='OR')
    ev <- subset(tcas,stratum=='EV')
    pv <- subset(tcas,stratum=='PV')
    sc <- subset(tcas,stratum=='SC')
    pt <- subset(tcas,stratum=='PT')
  
  for(i in 1:length(strata)){
    assign(paste0(strata[i],'01'),crop(preds01,get(tolower(strata[i]))))
    assign(paste0(strata[i],'01'),mask(get(paste0(strata[i],'01')),get(tolower(strata[i]))))
    assign(paste0(strata[i],'20'),crop(preds20,get(tolower(strata[i]))))
    assign(paste0(strata[i],'20'),mask(get(paste0(strata[i],'20')),get(tolower(strata[i]))))
  }
  
  tca.table <- data.frame(tca=strata)
  for(i in 1:length(strata)){
    tca.table$area.01[i] <- sum(!is.na(get(paste0(tca.table$tca[i],'01'))@data@values))
    tca.table$area.20[i] <- sum(!is.na(get(paste0(tca.table$tca[i],'20'))@data@values))
    tca.table$D.01[i] <- round(mean(get(paste0(tca.table$tca[i],'01'))@data@values,na.rm=TRUE),2)
    tca.table$D.20[i] <- round(mean(get(paste0(tca.table$tca[i],'20'))@data@values,na.rm=TRUE),2)
    tca.table$D.01.min[i] <- round(min(get(paste0(tca.table$tca[i],'01'))@data@values,na.rm=TRUE),2)
    tca.table$D.01.max[i] <- round(max(get(paste0(tca.table$tca[i],'01'))@data@values,na.rm=TRUE),2)
    tca.table$D.20.min[i] <- round(min(get(paste0(tca.table$tca[i],'20'))@data@values,na.rm=TRUE),2)
    tca.table$D.20.max[i] <- round(max(get(paste0(tca.table$tca[i],'20'))@data@values,na.rm=TRUE),2)
    tca.table$N.01[i] <- round(sum(get(paste0(tca.table$tca[i],'01'))@data@values,na.rm=TRUE))
    tca.table$N.20[i] <- round(sum(get(paste0(tca.table$tca[i],'20'))@data@values,na.rm=TRUE))
  }  
  tca.table
  

#-----------------------------------------------------------------------------------------------# 
# Plotting 2020 predicted densities
#-----------------------------------------------------------------------------------------------#     

  # Convert rasters to dataframes for ggplot
  
    preds20_df <- as.data.frame(preds20, xy = TRUE)

    imp19.list <- list(imp19.cd, imp19.em, imp19.nm, imp19.wm)
    imp19.binary <- do.call(merge, imp19.list)
    imp19.binary[imp19.binary >= is.threshold] <- 1
    imp19.binary[imp19.binary != 1] <- NA
    imp19_df <- as.data.frame(imp19.binary, xy = TRUE)    
  
  # Get range of values
    
    rng20 <- range(preds20_df$layer, na.rm = TRUE)
    midvalue20 <- mean(rng20)
    
  # Prep shapefiles for ggplot
    
    rus@data$id <- rownames(rus@data)
    rus_data <- fortify(rus, unit = "id")
    rus_df <- merge(rus_data, rus@data, by = "id")
    
    tcas_noUVR <- subset(tcas, recovery_u != "Upper Virgin River")
    tcas_noUVR@data$id <- rownames(tcas_noUVR@data)
    tcas_data <- fortify(tcas_noUVR, unit = "id")
    tcas_df <- merge(tcas_data, tcas_noUVR@data, by = "id")
    
    cstates <- subset(states, !STUSPS %in% c("MP", "AS", "PR", "VI", "HI" , "AK", "GU"))
    cstates@data$id <- rownames(cstates@data)
    cstates_data <- fortify(cstates, unit = "id")
    cstates_df <- merge(cstates_data, cstates@data, by = "id")    

    predarea@data$id <- rownames(predarea@data)
    predarea_data <- fortify(predarea, unit = "id")
    predarea_df <- merge(predarea_data, predarea@data, by = "id")

  # Plot
    
    cols <- c('lightsteelblue3','cornsilk','salmon3')
    options(scipen =1000000)
    
    g20 <- ggplot() + 
             geom_tile(data = preds20_df, aes(x = x, y = y, fill = layer)) + 
             scale_fill_gradient2(low = cols[1], mid = cols[2], high = cols[3], na.value = "white",
                                  midpoint = midvalue20, 
                                  guide = guide_colorbar(), 
                                  limits = rng20) +
             geom_path(data = cstates_df, aes(x = long, y = lat, group = group), 
                       color = "gray70", size = 0.8) +     
             geom_tile(data = dplyr::filter(imp19_df, !is.na(layer)), aes(x = x, y = y, fill = layer), 
                         fill = " gray30", show.legend = FALSE) +
             theme(axis.title = element_blank(),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   panel.background = element_blank(),
                   panel.border = element_rect(color = 'black', fill = NA),
                   legend.position = c(0.91, 0.45),
                   legend.title = element_blank(),
                   legend.key.height = unit(0.4, "in"),
                   legend.key.width = unit(0.2, "in"),
                   plot.margin = grid::unit(c(0.1,0,0,0), "in")) +
             geom_polygon(data = rus_df, aes(x = long, y = lat, group = group), 
                          fill = NA, color = 'black', size = 1) + 
             geom_polygon(data = tcas_df, aes(x = long, y = lat, group = group), 
                          fill = NA, color = 'black', size = 0.5) + 
             annotate("text", x = 780000, y = 3753000, label = "Density", hjust = 0, parse = TRUE, size = 11*0.8/.pt) + 
             annotate("text", x = 780000, y = 3740000, label = "Adults/km^2", hjust = 0, parse = TRUE, size = 11*0.8/.pt) + 
             coord_fixed(ratio = 1, xlim=c(336500, 825000), ylim=c(3631000, 4155000))
    g20
    
    ggsave("Figures/PredictedDensity_2020.jpg",
           device="jpeg",
           height = 6.8,
           width = 6.5,
           units = "in",
           dpi = 600)
    
