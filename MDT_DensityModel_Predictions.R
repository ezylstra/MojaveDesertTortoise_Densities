#===============================================================================================# 

 ## Authors: ER Zylstra, RJ Steidl

 ## Project: MDT line-distance sampling - estimating density (spatial variation and trends)

 ## Relevant files: 

  # 1. MDT_Import_Format_Data.R -- import all the files needed for analysis, format data, and create MDT_Data.Rdata
  # 2. MDT_DensityModel.R -- run model in STAN, summarize results, produce figures (trends)
  # 3. MDT_DensityModel.stan -- STAN model
  # 4. MDT_DensityModel_Predictions.R -- predict density of tortoises across their range (produce heat map)

#===============================================================================================# 

#-----------------------------------------------------------------------------------------------# 
# Set working directory, load packages
#-----------------------------------------------------------------------------------------------# 

  # Set working directory

    # setwd()

  # Load packages

    library(terra)
    library(tidyterra)
    library(ggplot2)
    library(gridExtra)
    library(dplyr)

#-----------------------------------------------------------------------------------------------# 
# Import data
#-----------------------------------------------------------------------------------------------# 
  
  # Import survey and associated covariate data 

    load('MDT_Data.Rdata')
    
  # Output from Stan model (and associated covariate data that went into the model)
    
    load('STANFit_MDT_Density_2001-2020.Rdata.Rdata')  
    posterior7500 <- as.matrix(out)

  # Covariate data    
    
    precip.s <- rast('Covariates/30yr_may_oct_prcp_mm.tif')            #Summer precip (May-Oct, mm), 30yr norms
    names(precip.s) <- 'precips'
    precip.w <- rast('Covariates/30yr_nov_apr_prcp_mm.tif')            #Winter precip (Nov-Apr, mm), 30yr norms
    names(precip.w) <- 'precipw'
    temp.max <- rast('Covariates/max_temp_warmest_month.tif')          #Maximum temperatures (degC), warmest month, 30yr norms
    names(temp.max) <- 'temp.max'
    north <- rast('Covariates/avg_northness.tif')                      #Aspect:northness 
    names(north) <- 'north'
    east <- rast('Covariates/avg_eastness.tif')                        #Aspect:eastness 
    names(east) <- 'east'
    rough <- rast('Covariates/average_surface_roughness_snapped.tif')  #Average surface roughness
    names(rough) <- 'rough'
    wash <- rast('Covariates/wash_proportion.tif')                     #Wash density (proportion of each pixel predicted to be "wash")
    names(wash) <- 'wash'
    bedrock <- rast('Covariates/bedrock_depth.tif')                    #Depth to bedrock (cm)
    names(bedrock) <- 'bedrock'
    veg.a <- rast('Covariates/annProx_snapped.tif')                    #Annual plant growth potential (MODIS, compare wet(2005)-dry(2002) years)
    names(veg.a) <- 'vega'
    road.a <- rast('Covariates/dist_road.tif')                         #Distance to any road (m)
    names(road.a) <- 'roada'

  # Shapefile: Recovery units
    
    rus <- vect('Covariates/Revised Recovery Units/2011RecoveryUnits.shp')
    rus <- project(rus, crs(east))

  # Shapefile: State boundaries
    
    states <- vect('Covariates/US_states_GIS/cb_2017_us_state_500k.shp')
    states <- project(states, crs(east))

  # Shapefile: TCAs
    
    tcas <- vect('Covariates/TCAs/All_Strata.shp')
    tcas <- project(tcas, crs(east))

  # Shapefile: Prediction Area
    # Includes any cells within 4 RUs that meet covariate criteria below
    # Excludes 1 sq km cells where one or more covariates values are NA
    # Excludes cells where any covariate is more than 10% outside the range of 
      # values observed at survey locations
    
    predarea <- vect('Covariates/PredictionArea/PredictionArea.shp')
    predarea <- project(predareafull, crs(east))

  # Shapefiles: major roads (from TIGER)
    
    roads.nv <- vect('Covariates/Roads_NV/tl_2021_32_prisecroads.shp')
    roads.nv <- project(roads.nv, crs(east))
    roads.ca <- vect('Covariates/Roads_CA/tl_2021_06_prisecroads.shp')
    roads.ca <- project(roads.ca, crs(east))
    roads.az <- vect('Covariates/Roads_AZ/tl_2021_04_prisecroads.shp')
    roads.az <- project(roads.az, crs(east))
    roads.ut <- vect('Covariates/Roads_UT/tl_2021_49_prisecroads.shp')
    roads.ut <- project(roads.ut, crs(east))
    highways.nv <- subset(roads.nv, roads.nv$MTFCC == 'S1100')
    highways.ca <- subset(roads.ca, roads.ca$MTFCC == 'S1100')
    highways.az <- subset(roads.az, roads.az$MTFCC == 'S1100')
    highways.ut <- subset(roads.ut, roads.ut$MTFCC == 'S1100')

  # Impervious surfaces
    
    # Raw values (%) for each 1 sq km cell in 2001
    imp01 <- rast('Covariates/ImpervSurfaces/impervious_raw_2001.tif')   
    
    # Raw values (%) for each 1 sq km cell in 2019
    imp19 <- rast('Covariates/ImpervSurfaces/impervious_raw_2019.tif')   

#-----------------------------------------------------------------------------------------------# 
# Extract survey locations
#-----------------------------------------------------------------------------------------------# 

  mids <- surveys[,c('segmentID','mid_easting','mid_northing')]
  midlocs <- vect(mids, geom = c("mid_easting", "mid_northing"), crs = crs(east))
 
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
    temp.mn <- mean(surveys$temp.max)
    rough.l.mn <- mean(surveys$rough.log)
    east.mn <- mean(surveys$aspect.e)
    north.mn <- mean(surveys$aspect.n)
    wash.l.mn <- mean(surveys$wash.log)
    bedrock.l.mn <- mean(surveys$bedrock.log)
    vega.mn <- mean(surveys$veg.a)
    roada.l.mn <- mean(surveys$road.a.log)
      
    precipw.sd <- sd(surveys$precip.w)
    precips.sd <- sd(surveys$precip.s)
    temp.sd <- sd(surveys$temp.max)
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
    temp.z <- (temp.max - temp.mn)/temp.sd
    rough.z <- (rough.log - rough.l.mn)/rough.l.sd
    east.z <- (east - east.mn)/east.sd
    north.z <- (north - north.mn)/north.sd
    wash.z <- (wash.log - wash.l.mn)/wash.l.sd
    bedrock.z <- (bedrock.log - bedrock.l.mn)/bedrock.l.sd
    vega.z <- (veg.a - vega.mn)/vega.sd
    roada.z <- (road.a.log - roada.l.mn)/roada.l.sd   
    
  # Create quadratic for temperature

    temp.z2 <- temp.z * temp.z
    names(temp.z2) <- "temp.max2"

#-----------------------------------------------------------------------------------------------#     
# Crop standardized covariate rasters to prediction areas within each recovery unit 
#-----------------------------------------------------------------------------------------------# 

  # Polygon for each RU
    
    cd <- subset(rus, rus$Unit_Name=='Colorado Desert')
    em <- subset(rus, rus$Unit_Name=='Eastern Mojave')
    nm <- subset(rus, rus$Unit_Name=='Northeastern Mojave')
    wm <- subset(rus, rus$Unit_Name=='Western Mojave')
  
  # Create shapefile with prediction area in each recovery unit
    
    pred.cd <- crop(predareafull, cd)
    pred.em <- crop(predareafull, em)
    pred.nm <- crop(predareafull, nm)
    pred.wm <- crop(predareafull, wm)

  # Crop (and mask) all covariate layers in CD
    
    precips.cd <- crop(precips.z,pred.cd)
    precips.cd <- mask(precips.cd,pred.cd)
    precipw.cd <- crop(precipw.z,pred.cd)
    precipw.cd <- mask(precipw.cd,pred.cd)
    temp.cd <- crop(temp.z,pred.cd)
    temp.cd <- mask(temp.cd,pred.cd)
    temp2.cd <- crop(temp.z2,pred.cd)
    temp2.cd <- mask(temp2.cd,pred.cd)    
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

  # Crop (and mask) all covariate layers in EM
    
    precips.em <- crop(precips.z,pred.em)
    precips.em <- mask(precips.em,pred.em)
    precipw.em <- crop(precipw.z,pred.em)
    precipw.em <- mask(precipw.em,pred.em)
    temp.em <- crop(temp.z,pred.em)
    temp.em <- mask(temp.em,pred.em)
    temp2.em <- crop(temp.z2,pred.em)
    temp2.em <- mask(temp2.em,pred.em)   
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

  # Crop (and mask) all covariate layers in NM
    
    precips.nm <- crop(precips.z,pred.nm)
    precips.nm <- mask(precips.nm,pred.nm)
    precipw.nm <- crop(precipw.z,pred.nm)
    precipw.nm <- mask(precipw.nm,pred.nm)
    temp.nm <- crop(temp.z,pred.nm)
    temp.nm <- mask(temp.nm,pred.nm)
    temp2.nm <- crop(temp.z2,pred.nm)
    temp2.nm <- mask(temp2.nm,pred.nm)   
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

  # Crop (and mask) all covariate layers in WM
    
    precips.wm <- crop(precips.z,pred.wm)
    precips.wm <- mask(precips.wm,pred.wm)
    precipw.wm <- crop(precipw.z,pred.wm)
    precipw.wm <- mask(precipw.wm,pred.wm)
    temp.wm <- crop(temp.z,pred.wm)
    temp.wm <- mask(temp.wm,pred.wm)
    temp2.wm <- crop(temp.z2,pred.wm)
    temp2.wm <- mask(temp2.wm,pred.wm)   
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

#-----------------------------------------------------------------------------------------------# 
# Create impervious surface layers
#-----------------------------------------------------------------------------------------------#

  is.threshold <- 0.4

  # Create binary layers (where 1/TRUE indicates impervious values > threshold)
    
    impB01 <- imp01 > is.threshold
    impB19 <- imp19 > is.threshold
  
  # Create layers where 1 indicates cells that have impervious values < threshold
  # and all other cells are NA
    
    nimp01 <- impB01
    nimp01[nimp01 == 1] <- NA
    nimp01[nimp01 == 0] <- 1
    nimp19 <- impB19
    nimp19[nimp19 == 1] <- NA
    nimp19[nimp19 == 0] <- 1
  
  # Crop (and mask) impervious surface layers in CD
    
    nimp01.cd <- crop(nimp01, pred.cd)
    nimp01.cd <- mask(nimp01.cd, pred.cd)
    nimp19.cd <- crop(nimp19, pred.cd)
    nimp19.cd <- mask(nimp19.cd, pred.cd)

  # Crop (and mask) impervious surface layers in EM
    
    nimp01.em <- crop(nimp01, pred.em)
    nimp01.em <- mask(nimp01.em, pred.em)
    nimp19.em <- crop(nimp19, pred.em)
    nimp19.em <- mask(nimp19.em, pred.em)
  
  # Crop (and mask) impervious surface layers in NM
    
    nimp01.nm <- crop(nimp01, pred.nm)
    nimp01.nm <- mask(nimp01.nm, pred.nm)
    nimp19.nm <- crop(nimp19, pred.nm)
    nimp19.nm <- mask(nimp19.nm, pred.nm)
  
  # Crop (and mask) impervious surface layers in WM
    
    nimp01.wm <- crop(nimp01, pred.wm)
    nimp01.wm <- mask(nimp01.wm, pred.wm)
    nimp19.wm <- crop(nimp19, pred.wm)
    nimp19.wm <- mask(nimp19.wm, pred.wm)

#-----------------------------------------------------------------------------------------------# 
# Organize covariates for prediction in each grid cell
#-----------------------------------------------------------------------------------------------#  
    
  # Putting covariates in the same order they appear in the model (cov_lam)
    
    covs.cd <- rast(list(east.cd,north.cd,bedrock.cd,precips.cd,precipw.cd,
                         roada.cd,rough.cd,temp.cd,temp2.cd,vega.cd,wash.cd))
    covs.em <- rast(list(east.em,north.em,bedrock.em,precips.em,precipw.em,
                         roada.em,rough.em,temp.em,temp2.em,vega.em,wash.em))
    covs.nm <- rast(list(east.nm,north.nm,bedrock.nm,precips.nm,precipw.nm,
                         roada.nm,rough.nm,temp.nm,temp2.nm,vega.nm,wash.nm))
    covs.wm <- rast(list(east.wm,north.wm,bedrock.wm,precips.wm,precipw.wm,
                         roada.wm,rough.wm,temp.wm,temp2.wm,vega.wm,wash.wm))

  # Multiply by nimp layers, so any cells with impervious surface values above 
  # the threshold become NA
    
    covs01.cd <- covs.cd * nimp01.cd
    covs19.cd <- covs.cd * nimp19.cd
    covs01.em <- covs.em * nimp01.em
    covs19.em <- covs.em * nimp19.em
    covs01.nm <- covs.nm * nimp01.nm
    covs19.nm <- covs.nm * nimp19.nm
    covs01.wm <- covs.wm * nimp01.wm
    covs19.wm <- covs.wm * nimp19.wm
    
  # Convert each raster to a dataframe:
  # Note: each column is a layer, rows are in order of raster cells
  
    covm01.cd <- as.data.frame(covs01.cd, cells = TRUE, na.rm = TRUE)
    covm19.cd <- as.data.frame(covs19.cd, cells = TRUE, na.rm = TRUE)
    covm01.em <- as.data.frame(covs01.em, cells = TRUE, na.rm = TRUE)
    covm19.em <- as.data.frame(covs19.em, cells = TRUE, na.rm = TRUE)
    covm01.nm <- as.data.frame(covs01.nm, cells = TRUE, na.rm = TRUE)
    covm19.nm <- as.data.frame(covs19.nm, cells = TRUE, na.rm = TRUE)
    covm01.wm <- as.data.frame(covs01.wm, cells = TRUE, na.rm = TRUE)
    covm19.wm <- as.data.frame(covs19.wm, cells = TRUE, na.rm = TRUE)

#-----------------------------------------------------------------------------------------------# 
# Calculating predicted densities for each raster cell in 2001 (after removing impervious surfaces)
#-----------------------------------------------------------------------------------------------#     

  # Note: Not including random yearly effects around range-wide trend, which are more likely to 
  # capture differences in areas sampled each year than true density changes

  # Extracting posterior samples for regression coefficients in abundance model
    # Note: order of covariates in rasters and in model parameters (cov_lam) must be the same
    # east, north, bedrock, precips, precipw, road, rough, temp, temp2, veg, wash

    iter <- floor(seq(1, 7500, length = 1000))
    
    # Regression coefficients in the abundance model (on log scale)
      
      betas <- as.matrix(posterior7500[iter,grep('beta_lam',colnames(posterior7500))])
    
    # Recovery unit intercepts (density on real scale in year 2000)
      
      mus <- as.matrix(posterior7500[iter,grep('mu_recov',colnames(posterior7500))])
    
    # Regression coefficient for yearly trends in the abundance model (on log scale)
      
      trends <- as.matrix(posterior7500[iter,grep('trend',colnames(posterior7500))])

  # Calculating predicted densities in 2001
    # log(pred.yr1) = log(mus[RU]) + 1 %*% trend[RU] + covm01.ru %*% betas
    # order of RUs: CD, EM, NM, WM
    
    # CD
      
      cd.mu <- as.matrix(log(mus[,1]))
      cd.trend <- as.matrix(trends[,1])
      cd.predl01 <- matrix(1,nrow=nrow(covm01.cd),ncol=1) %*% t(cd.mu) + 
                    matrix(1,nrow=nrow(covm01.cd),ncol=1) %*% t(cd.trend) +
                    as.matrix(covm01.cd[,-1]) %*% t(betas)
      cd.pred01 <- exp(cd.predl01)
      # Median expected density in each cell
      cd.predm01 <- apply(cd.pred01,1,median) 
      # Creating raster with median density values
      cd.preds01 <- rast(covs01.cd[[1]])
      cd.preds01[covm01.cd[,1]] <- cd.predm01

    # EM
      
      em.mu <- as.matrix(log(mus[,2]))
      em.trend <- as.matrix(trends[,2]) 
      em.predl01 <- matrix(1,nrow=nrow(covm01.em),ncol=1) %*% t(em.mu) + 
                    matrix(1,nrow=nrow(covm01.em),ncol=1) %*% t(em.trend) +
                    as.matrix(covm01.em[,-1]) %*% t(betas)
      em.pred01 <- exp(em.predl01)
      # Median expected density in each cell
      em.predm01 <- apply(em.pred01,1,median) 
      # Creating raster with median density values
      em.preds01 <- rast(covs01.em[[1]])
      em.preds01[covm01.em[,1]] <- em.predm01

    # NM
      
      nm.mu <- as.matrix(log(mus[,3]))
      nm.trend <- as.matrix(trends[,3])      
      nm.predl01 <- matrix(1,nrow=nrow(covm01.nm),ncol=1) %*% t(nm.mu) + 
                    matrix(1,nrow=nrow(covm01.nm),ncol=1) %*% t(nm.trend) +
                    as.matrix(covm01.nm[,-1]) %*% t(betas)
      nm.pred01 <- exp(nm.predl01)
      # Median expected density in each cell
      nm.predm01 <- apply(nm.pred01,1,median) 
      # Creating raster with median density values
      nm.preds01 <- rast(covs01.nm[[1]])
      nm.preds01[covm01.nm[,1]] <- nm.predm01

    # WM
      
      wm.mu <- as.matrix(log(mus[,4]))
      wm.trend <- as.matrix(trends[,4])      
      wm.predl01 <- matrix(1,nrow=nrow(covm01.wm),ncol=1) %*% t(wm.mu) + 
                    matrix(1,nrow=nrow(covm01.wm),ncol=1) %*% t(wm.trend) +
                    as.matrix(covm01.wm[,-1]) %*% t(betas)
      wm.pred01 <- exp(wm.predl01)
      # Median expected density in each cell
      wm.predm01 <- apply(wm.pred01,1,median) 
      # Creating raster with median density values
      wm.preds01 <- rast(covs01.wm[[1]])
      wm.preds01[covm01.wm[,1]] <- wm.predm01

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
                    as.matrix(covm19.cd[,-1]) %*% t(betas)
      cd.pred20 <- exp(cd.predl20)
      cd.predm20 <- apply(cd.pred20,1,median)
      cd.preds20 <- rast(covs19.cd[[1]])
      cd.preds20[covm19.cd[,1]] <- cd.predm20

    # EM
      
      em.predl20 <- matrix(1,nrow=nrow(covm19.em),ncol=1) %*% t(em.mu) + 
                    matrix(20,nrow=nrow(covm19.em),ncol=1) %*% t(em.trend) +
                    as.matrix(covm19.em[,-1]) %*% t(betas)
      em.pred20 <- exp(em.predl20)
      em.predm20 <- apply(em.pred20,1,median)
      em.preds20 <- rast(covs19.em[[1]])
      em.preds20[covm19.em[,1]] <- em.predm20

    # NM
      
      nm.predl20 <- matrix(1,nrow=nrow(covm19.nm),ncol=1) %*% t(nm.mu) + 
                    matrix(20,nrow=nrow(covm19.nm),ncol=1) %*% t(nm.trend) +
                    as.matrix(covm19.nm[,-1]) %*% t(betas)
      nm.pred20 <- exp(nm.predl20)
      nm.predm20 <- apply(nm.pred20,1,median)
      nm.preds20 <- rast(covs19.nm[[1]])
      nm.preds20[covm19.nm[,1]] <- nm.predm20

    #WM
      
      wm.predl20 <- matrix(1,nrow=nrow(covm19.wm),ncol=1) %*% t(wm.mu) + 
                    matrix(20,nrow=nrow(covm19.wm),ncol=1) %*% t(wm.trend) +
                    as.matrix(covm19.wm[,-1]) %*% t(betas)
      wm.pred20 <- exp(wm.predl20)
      wm.predm20 <- apply(wm.pred20,1,median)
      wm.preds20 <- rast(covs19.wm[[1]])
      wm.preds20[covm19.wm[,1]] <- wm.predm20

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
    area.cd01 <- nrow(cd.pred01)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.cd01 <- apply(cd.pred01,2,sum)
    A.cd01 <- round(median(A1000.cd01))
    A.cd01.sd <- round(sd(A1000.cd01))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.cd01 <- apply(cd.pred01,2,mean)
    D.cd01 <- round(median(D1000.cd01),2)
    D.cd01.sd <- round(sd(D1000.cd01),2)

  # CD - 2020
    
    # Area
    area.cd20 <- nrow(cd.pred20)
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
    area.em01 <- nrow(em.pred01)
    # Abundance (for each iteration, sum densities across cells; calculate median across 1000 sums)
    A1000.em01 <- apply(em.pred01,2,sum)
    A.em01 <- round(median(A1000.em01))
    A.em01.sd <- round(sd(A1000.em01))    
    # Density (for each iteration, average densities across cells; calculate median across 1000 means)
    D1000.em01 <- apply(em.pred01,2,mean)
    D.em01 <- round(median(D1000.em01),2)
    D.em01.sd <- round(sd(D1000.em01),2)

  # EM - 2020
    
    # Area
    area.em20 <- nrow(em.pred20)
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
    area.nm01 <- nrow(nm.pred01)
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
    area.nm20 <- nrow(nm.pred20)
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
    area.wm01 <- nrow(wm.pred01)
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
    area.wm20 <- nrow(wm.pred20)
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
  
    strata <- as.data.frame(tcas)$stratum
    strata <- strata[strata != "RC"]
    ag <- subset(tcas, tcas$stratum=='AG')
    jt <- subset(tcas, tcas$stratum=='JT')
    bd <- subset(tcas, tcas$stratum=='BD')
    ck <- subset(tcas, tcas$stratum=='CK')
    cm <- subset(tcas, tcas$stratum=='CM')
    cs <- subset(tcas, tcas$stratum=='CS')
    fe <- subset(tcas, tcas$stratum=='FE')
    fk <- subset(tcas, tcas$stratum=='FK')
    gb <- subset(tcas, tcas$stratum=='GB')
    iv <- subset(tcas, tcas$stratum=='IV')
    mm <- subset(tcas, tcas$stratum=='MM')
    or <- subset(tcas, tcas$stratum=='OR')
    ev <- subset(tcas, tcas$stratum=='EV')
    pv <- subset(tcas, tcas$stratum=='PV')
    sc <- subset(tcas, tcas$stratum=='SC')
    pt <- subset(tcas, tcas$stratum=='PT')
  
  for(i in 1:length(strata)){
    assign(paste0(strata[i],'01'),crop(preds01,get(tolower(strata[i]))))
    assign(paste0(strata[i],'01'),mask(get(paste0(strata[i],'01')),get(tolower(strata[i]))))
    assign(paste0(strata[i],'20'),crop(preds20,get(tolower(strata[i]))))
    assign(paste0(strata[i],'20'),mask(get(paste0(strata[i],'20')),get(tolower(strata[i]))))
  }
  
  tca.table <- data.frame(tca = strata)
  for(i in 1:length(strata)){
    tca.table$area.01[i] <- unlist(global(get(paste0(tca.table$tca[i],'01')), fun = "notNA"))
    tca.table$area.20[i] <- unlist(global(get(paste0(tca.table$tca[i],'20')), fun = "notNA"))
    tca.table$D.01[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'01')), fun = "mean", na.rm=TRUE),2))
    tca.table$D.20[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'20')), fun = "mean", na.rm=TRUE),2))
    tca.table$D.01.min[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'01')), fun = "min", na.rm=TRUE),2))
    tca.table$D.01.max[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'01')), fun = "max", na.rm=TRUE),2))
    tca.table$D.20.min[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'20')), fun = "min", na.rm=TRUE),2))
    tca.table$D.20.max[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'20')), fun = "max", na.rm=TRUE),2))
    tca.table$N.01[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'01')), fun = "sum", na.rm=TRUE)))
    tca.table$N.20[i] <- unlist(round(global(get(paste0(tca.table$tca[i],'20')), fun = "sum", na.rm=TRUE)))
  }
  tca.table

#-----------------------------------------------------------------------------------------------# 
# Plotting predicted densities in 2020
#-----------------------------------------------------------------------------------------------#     

# Get range of values just for 2020
  
  rng20 <- as.numeric(global(preds20, fun = "range", na.rm = TRUE))
  midvalue20 <- mean(rng20)

# Format shapefiles 
  
  tcas_noUVR <- subset(tcas, tcas$stratum != "RC")
  cstates <- subset(states, !states$STUSPS %in% c("MP", "AS", "PR", "VI", "HI" , "AK", "GU"))
  azI15 <- subset(roads.az, roads.az$FULLNAME == 'I- 15')

# Save layers in a different projection for study area map inset
  
  newcrs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  cstates_albers <- terra::project(cstates, newcrs)
  rus_albers <- terra::project(rus, newcrs)

# Prep objects for ggplot
  
  cstates <- ggplot2::fortify(cstates)
  rus <- ggplot2::fortify(rus)
  tcas_noUVR <- ggplot2::fortify(tcas_noUVR)
  preds20 <- ggplot2::fortify(preds20)
  
  impB19_plot <- impB19
  names(impB19_plot) <- "imp"
  impB19_plot[impB19_plot != 1] <- NA
  impB19_plot <- crop(impB19_plot, predareafull)
  impB19_plot <- mask(impB19_plot, predareafull)
  impB19_plot <- ggplot2::fortify(impB19_plot)
  

# Create plot

  cols <- c('lightsteelblue3','cornsilk','salmon3')  
    
  g20 <- ggplot() +
    geom_tile(data = preds20, aes(x = x, y = y, fill = east))  +
    scale_fill_gradient2(low = cols[1], mid = cols[2], high = cols[3], na.value = "white",
                         midpoint = midvalue20, 
                         guide = guide_colorbar(), 
                         limits = rng20) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.background = element_blank(),
          panel.border = element_rect(color = 'black', fill = NA),
          legend.position = c(0.91, 0.45),
          legend.title = element_blank(),
          legend.key.height = unit(0.4, "in"),
          legend.key.width = unit(0.2, "in"),
          plot.margin = grid::unit(c(0.1,0,0,0), "in")) +
    geom_tile(data = impB19_plot[!is.na(impB19_plot$imp),], 
              aes(x = x, y = y, fill = imp), 
              fill = "gray30", show.legend = FALSE) +
    geom_sf(data = cstates, fill = NA, color = "gray70", linewidth = 0.5) +
    geom_sf(data = rus, fill = NA, color = 'black', linewidth = 0.7) +
    geom_sf(data = tcas_noUVR, fill = NA, color = 'black', linewidth = 0.4) +
    annotate("text", x = 780000, y = 3753000, label = "Density", hjust = 0, parse = TRUE, size = 11*0.8/.pt) + 
    annotate("text", x = 780000, y = 3740000, label = "Adults/km^2", hjust = 0, parse = TRUE, size = 11*0.8/.pt) +          
    coord_sf(xlim=c(336500, 825000), ylim=c(3631000, 4155000))
  
  g20
  
  ggsave("Figures/PredictedDensity_2020.pdf",
         device="pdf",
         height = 6.8,
         width = 6.5,
         units = "in")

#-----------------------------------------------------------------------------------------------# 
# Study area plot
#-----------------------------------------------------------------------------------------------#   

  plotlocs <- unique(g0obs[,c("site","east","north")])
  plotlocs <- vect(plotlocs, geom = c("east", "north"), crs = crs(east))

# Prep objects for ggplot

  predareafull <- ggplot2::fortify(predareafull)
  midlocs <- ggplot2::fortify(midlocs)
  plotlocs <- ggplot2::fortify(plotlocs)
  highways.ca <- ggplot2::fortify(highways.ca)
  highways.ut <- ggplot2::fortify(highways.ut)
  highways.nv <- ggplot2::fortify(highways.nv)
  highways.az <- ggplot2::fortify(highways.az)
  cstates_albers <- ggplot2::fortify(cstates_albers)
  rus_albers <- ggplot2::fortify(rus_albers)

# Create plot
  
  gSA <- ggplot() +
    geom_sf(data = predareafull, fill = "gray85", color = NA) + 
    geom_sf(data = midlocs, color = "steelblue3", size = 0.3) + 
    geom_sf(data = highways.nv, color = "gray50", linewidth = 0.05) +
    geom_sf(data = highways.ca, color = "gray50", linewidth = 0.05) +
    geom_sf(data = highways.ut, color = "gray50", linewidth = 0.05) +
    geom_sf(data = highways.az, color = "gray50", linewidth = 0.05) +
    geom_sf(data = cstates, fill = NA, color = "gray50", linewidth = 0.6) +
    geom_sf(data = rus, fill = NA, color = "black", linewidth = 0.7) +
    geom_sf(data = tcas_noUVR, fill = NA, col = "black", linewidth = 0.5) +
    geom_sf(data = plotlocs, fill = "yellow", color = "black", pch = 23, size = 2) +
    annotate("text", x = 470000, y = 3970000, label = "WM", hjust = 0.5, vjust = 0.5, size = 12/.pt, fontface = 2) +
    annotate("text", x = 580000, y = 4080000, label = "EM", hjust = 0.5, vjust = 0.5, size = 12/.pt, fontface = 2) +
    annotate("text", x = 700000, y = 4138000, label = "NM", hjust = 0.5, vjust = 0.5, size = 12/.pt, fontface = 2) +
    annotate("text", x = 700000, y = 3752000, label = "CD", hjust = 0.5, vjust = 0.5, size = 12/.pt, fontface = 2) +
    annotate("text", x = 805000, y = 4120000, label = "UVR", hjust = 0.5, vjust = 0.5, size = 12/.pt, fontface = 2) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.background = element_blank(),
          panel.border = element_rect(color = 'black', fill = NA),
          plot.margin = grid::unit(c(0.1,0,0,0), "in")) + 
    coord_sf(xlim=c(336500, 825000), ylim=c(3631000, 4155000))
  
  gSA
  
  ginset <- ggplot(cstates_albers) + 
    geom_sf(fill = NA, color = "gray30", size = 0.2) + 
    geom_sf(data = rus_albers, fill = "black") + 
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = 'black', fill = NA))
  
  gSA + annotation_custom(grob=ggplotGrob(ginset), 
                          xmin = 310000, xmax = 510000, ymin = 3600000, ymax = 3760000) 
  
  ggsave("Figures/StudyArea_wInset.pdf",
         device="pdf",
         height = 6.8,
         width = 6.5,
         units = "in")
