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

  # rm(list=ls())
  
#-----------------------------------------------------------------------------------------------# 
# Import data
#-----------------------------------------------------------------------------------------------# 
  
  #Distance-sampling data
    
    #Survey data (e.g., segment names, lengths)
    allSurveys <- read.csv('TortoiseData/Tort_Surveys.csv',header=TRUE,na.strings=c("NA",""),strip.white=TRUE)   
    #Detection data for live tortoises
    live <- read.csv('TortoiseData/Tort_Obs_Live.csv',header=TRUE,na.strings=c("NA",""),strip.white=TRUE) 
    #Recovery units associated with each stratum
    ru <- read.csv('TortoiseData/RUs.csv',header=TRUE,strip.white=TRUE)
    
  #Covariate data
    
    precip.s <- raster('Covariates/30yr_may_oct_prcp_mm.tif')            #Summer precip (May-Oct, mm), 30yr norms
    precip.w <- raster('Covariates/30yr_nov_apr_prcp_mm.tif')            #Winter precip (Nov-Apr, mm), 30yr norms
    temp.max <- raster('Covariates/max_temp_warmest_month.tif')          #Maximum temperatures (degC), warmest month, 30yr norms
    elev <- raster('Covariates/avg_elevation.tif')                       #Elevation (m)
    slope <- raster('Covariates/avg_slope.tif')                          #Slope (degrees)
    north <- raster('Covariates/avg_northness.tif')                      #Aspect:northness 
    east <- raster('Covariates/avg_eastness.tif')                        #Aspect:eastness 
    rough <- raster('Covariates/average_surface_roughness_snapped.tif')  #Average surface roughness
    wash <- raster('Covariates/wash_proportion.tif')                     #Wash density (proportion of each pixel predicted to be "wash")
    bedrock <- raster('Covariates/bedrock_depth.tif')                    #Depth to bedrock (cm)
    soilbulk <- raster('Covariates/soil_bulk_density.tif')               #Soil bulk density
    prock <- raster('Covariates/pct_rocks_gt254mm.tif')                  #Percent rock >254 mm
    veg.p <- raster('Covariates/perennial_vegetation.tif')               #Perennial plant cover: NDVI values in 2019 (15 Jun - 31 Aug), very dry year 
    veg.a <- raster('Covariates/annProx_snapped.tif')                    #Annual plant growth potential (MODIS, compare wet(2005)-dry(2002) years)
    road.m <- raster('Covariates/dist_minor_road.tif')                   #Distance to minor road (m)
    road.a <- raster('Covariates/dist_road.tif')                         #Distance to any road (m)
    #Winter precipitation data (Oct-Mar preceding each survey, mm)
    ppt.surveys <- read.csv('Covariates/survey_start_points_with_precip_DAYMET.csv',header=TRUE,strip.white=TRUE)
    ppt.octmar <- raster('Covariates/30yr_oct_mar_prcp_mm.tif')          #Winter precip (Oct-Mar, mm), 30yr norms

  #Telemetry data (for estimating availability/g0)
    
    #Telemetry data, with binary "visible" variable (excluding RM site)
    g0obs <- read.csv('TortoiseData/g0_Obs.csv',header=TRUE,na.strings=c("NA",""),strip.white=TRUE)        
    #Telemetry data from RM (teaching) site, 2012-2013 only
    g0obs.rm <- read.csv('TortoiseData/g0_Obs_RMonly.csv',header=FALSE,na.strings=c("NA",""),strip.white=TRUE,col.names=names(g0obs))                                                    
    #Name, abbreviation, and location of each telemetry site
    g0sites <- read.csv('TortoiseData/g0_Sites.csv',header=TRUE,na.strings=c("NA",""),strip.white=TRUE)
    
    #Winter precipitation data (Oct-Mar preceding each telemetry season, mm)
    ppt.g0 <-read.csv('Covariates/telemetry_sites_yearly_precip_DAYMET.csv',header=TRUE,strip.white=T)

#-----------------------------------------------------------------------------------------------# 
# Remove data from Upper Virgin River (UVR) recovery unit
#-----------------------------------------------------------------------------------------------#  
    
  #Add recovery unit to survey data, live observations
    
    allSurveys$recovID <- ru$recovery_unit[match(allSurveys$stratum,ru$stratum)]
    live$recovID       <- allSurveys$recovID[match(live$transectID,allSurveys$transectID)]

  #Remove distance-sampling data from UVR
    
    allSurveys <- allSurveys[allSurveys$recovID!='Upper Virgin River',]
    live       <- live[which(live$recovID!='Upper Virgin River'),]

  #Remove telemetry data from UVR (site = RC)
  
    g0obs <- g0obs[g0obs$G0_site!='RC',]
      
#-----------------------------------------------------------------------------------------------# 
# Establish start and end years
#-----------------------------------------------------------------------------------------------#  
  
  minYear <- 2001
  maxYear <- 2018

#-----------------------------------------------------------------------------------------------# 
# Format telemetry data
#-----------------------------------------------------------------------------------------------#  
    
  #Merge RM data with the rest of telemetry data
  
    g0obs <- rbind(g0obs,g0obs.rm)
  
  #Format date, year, and day-of-year
    
    g0obs$obsdate <- as.Date(g0obs$date_,format="%m/%d/%Y")  
    #Recreate year from date to avoid data entry errors (and name it "season" to match precipitation data)
    g0obs$season  <- as.numeric(format(g0obs$obsdate,'%Y'))   
    g0obs$doy <- as.numeric(format(g0obs$obsdate,'%j'))

  #Create unique tortoise IDs because tort_num is not unique across sites
  
    g0obs$tortid <- paste(g0obs$G0_site,g0obs$tort_num,sep='_')  
    
  #Convert response (visible) to numeric values
  
    # count(g0obs$visible)
    g0obs <- g0obs[g0obs$visible!='Unk',] 
    g0obs$visible[g0obs$visible %in% c('N','No')]  <- 0
    g0obs$visible[g0obs$visible %in% c('Y','Yes')] <- 1
    g0obs$visible <- as.numeric(g0obs$visible)
  
  #Make sure that all sites in g0obs are in g0sites
  
    #Should return "character(0)" if all sites in g0obs are in g0sites
    # unique(g0obs$G0_site[!g0obs$G0_site %in% g0sites$G0_site]) 
    
  #Add site UTMs and Oct-Mar precipitation normals to g0obs

    g0sites <- g0sites[,c('G0_site','Easting','Northing')]
    g0locs <- SpatialPointsDataFrame(coords=g0sites[,2:3],proj4string=crs(ppt.octmar),data=g0sites)
    g0sites$ppt.norm <- extract(ppt.octmar,g0locs)
    g0obs <- join(g0obs,g0sites,by='G0_site',type='left')
    
  #Add yearly precipitation data to g0obs
    
    #Check that both sets of data include the same sites
    # all.equal(sort(ppt.g0$g0site),sort(unique(g0obs$G0_site)))
    #Clean up yearly precipitation data and put in long form
    ppt.g0 <- ppt.g0[,c(1,5:ncol(ppt.g0))]
    names(ppt.g0) <- c('G0_site',paste('ppt',2001:2020,sep='.'))
    ppt.g0.long <- melt(ppt.g0,id.vars='G0_site',variable.name='ppt.yr',value.name='ppt')
    ppt.g0.long$season <- as.numeric(substr(ppt.g0.long$ppt.yr,5,8))
    #Add to g0obs
    g0obs <- join(g0obs,ppt.g0.long[,c('G0_site','ppt','season')],by=c('G0_site','season'),type='left')
    
  #Calculate annual % of 30-yr Oct-Mar precipitation norms 
    
    g0obs$pptOM.perc <- g0obs$ppt/g0obs$ppt.norm*100
    # summary(g0obs$pptOM.perc) #ranges from 7-315%

  #Keep only the columns we'll need and rename
    
    g0obs <- g0obs[,c('tortid','G0_site','Easting','Northing','obsdate','season','doy','visible','ppt','ppt.norm','pptOM.perc')]
    names(g0obs)[c(2:4)] <- c('site','east','north')  
    
  #Subset g0 data by year
    
    g0obs <- g0obs[g0obs$season %in% minYear:maxYear,]

#-----------------------------------------------------------------------------------------------# 
# Format observations of live tortoises (each row in dataframe represents one tortoise)
#-----------------------------------------------------------------------------------------------# 

  #Round distances to 1 decimal place 
    
    live$perp_distance_m <- round(live$perp_distance_m,1) 

  #Fix indicator, so mcl_greater_180='yes' when mcl measurement >= 180 mm:
    
    # count(live$mcl_greater_180)
    # summary(live$mcl_mm)
    # summary(live$mcl_mm[live$mcl_greater_180=='unknown'])  # All entries with mcl_greater_180=='unknown' have mcl_mm=NA
    # ddply(live[live$mcl_greater_180!='unknown',],.(mcl_greater_180),
    #       summarize,min=min(mcl_mm,na.rm=T),max=max(mcl_mm,na.rm=T),nNA=sum(is.na(mcl_mm)))
    live$mcl_greater_180[which(live$mcl_mm>=180)] <- 'yes'
    
  #Remove unnecessary variables (keeping mcl variables and sex so that we could subset data by either) and rename columns
      
    live <- live[,names(live) %in% c('segmentID','year_','perp_distance_m','sex','mcl_mm','mcl_greater_180'),]
    names(live)[c(1,3)] <- c('season','dist')
    
  #Subset observations by year
    
    live <- live[live$season %in% minYear:maxYear,]
    
  #Subset observations to exclude tortoises with MCL<180 (or exclude if we don't know whether MCL>=180)

    live <- live[live$mcl_greater_180=='yes',]  
    
  #Truncate observations by detection distance
    #If truncating, set Trunc = W = truncation distance (Eguchi and Gerrodette 2009)
    #If not truncating, set W = maxDist (max. perpendicular distance). Note: need to change the likelihood in STAN model if we want to do this!
    
    Trunc <- W <- 20   
    # nrow(live) #4962
    live <- live[live$dist<=Trunc,]
    # nrow(live) #4746 (truncating 4.4% of observations)
    
    liveyr <- ddply(live,.(season),summarize,ndetects=length(dist))
    summary(liveyr)
    
    live <- live[,c('segmentID','dist')]
    
#-----------------------------------------------------------------------------------------------# 
# Format survey data (each row in dataframe represents one segment)
#-----------------------------------------------------------------------------------------------#  
 
  #Remove unnecessary columns (including covariate data, which we'll get from raster layers)
    include <- c('year','stratum','transectID','segmentID','date_','seg_length','group','start_lead',
                 'start_east','start_nort','end_eastin','end_northi','recovID')
    allSurveys <- allSurveys[,names(allSurveys) %in% include]
    names(allSurveys) <- c('season','stratum','transectID','segmentID','date','seg_length','group',
                           'start_lead','start_e','start_n','end_e','end_n','recovID')

  #Format date, calculate day-of-year
    
    allSurveys$sdate <- as.Date(allSurveys$date,format='%m/%d/%Y')
    allSurveys$doy <- as.numeric(format(allSurveys$sdate,'%j'))
    allSurveys$date <- NULL
    
  #Subset survey data by year
    
    allSurveys <- allSurveys[allSurveys$season %in% minYear:maxYear,]
  
#-----------------------------------------------------------------------------------------------# 
# Add Oct-Mar precipitation data (ppt.surveys) to the survey data
#-----------------------------------------------------------------------------------------------#    

  #First, check that all segments in allSurveys are in the ppt.surveys dataframe
    
    #Should return "character(0)" if all segments in allSurveys are in ppt.surveys
    # unique(allSurveys$segmentID[!allSurveys$segmentID %in% ppt.surveys$segmentID])

  #Add annual precipitation data to survey data
    
    ppt.surveys <- ppt.surveys[,c('prev_6mo_precip','segmentID')]
    names(ppt.surveys)[1] <- 'ppt'
    allSurveys$ppt <- ppt.surveys$ppt[match(allSurveys$segmentID,ppt.surveys$segmentID)]
    
  #Add Oct-Mar precipitation normals to survey data
    
    starts <- allSurveys[,c('segmentID','start_e','start_n')]
    startlocs <- SpatialPointsDataFrame(coords=starts[,2:3],proj4string=crs(elev),data=starts)
    allSurveys$ppt.norm <- extract(ppt.octmar,startlocs)
    
  #Calculate annual % of 30-yr Oct-Mar precipitation norms 
    
    allSurveys$pptOM.perc <- allSurveys$ppt/allSurveys$ppt.norm*100
    # summary(allSurveys$pptOM.perc) #ranges from 5-334%  

#-----------------------------------------------------------------------------------------------# 
# Attach spatial covariate data to allSurveys
#-----------------------------------------------------------------------------------------------#   
    
  #The wash layer is missing some values in the north-central part of the study area
  #Imputing the mean value
  
    wash.impute <- wash
    wash.impute[is.na(wash.impute)] <- mean(values(wash),na.rm=TRUE)
    
  #Extract values for each segment based on start location
    
    starts$precip.s <- extract(precip.s,startlocs)
    starts$precip.w <- extract(precip.w,startlocs)
    starts$temp.max <- extract(temp.max,startlocs)
    starts$elev <- extract(elev,startlocs)
    starts$slope <- extract(slope,startlocs)
    starts$aspect.n <- extract(north,startlocs)
    starts$aspect.e <- extract(east,startlocs)
    starts$rough <- extract(rough,startlocs)
    starts$wash <- extract(wash.impute,startlocs)
    starts$bedrock <- extract(bedrock,startlocs)
    starts$soilbulk <- extract(soilbulk,startlocs)
    starts$prock <- extract(prock,startlocs)
    starts$veg.p <- extract(veg.p,startlocs)
    starts$veg.a <- extract(veg.a,startlocs)
    starts$road.m <- extract(road.m,startlocs)
    starts$road.a <- extract(road.a,startlocs)
    
    #Check that there are no missing values:
    summary(starts)

  #Merge covariate data with survey data

    allSurveys <- join(allSurveys,starts[,c(1,4:ncol(starts))],by='segmentID',type='left')
    
#-----------------------------------------------------------------------------------------------# 
# Merge survey and observation data
#-----------------------------------------------------------------------------------------------# 

  #Count the number of observations per segment
    
    nDetects <- as.data.frame(table(live$segmentID)) 
    names(nDetects) <- c('segmentID','nDetects')
    
  #Reshape detection distances from long to wide (one observation per row to all observations for a segment on one row)
    
    live <- live[with(live,order(segmentID,dist)),]      
    #Add observation number (1:nDetects) for each detection on each segment
    live$obs <- sequence(rle(as.character(live$segmentID))$lengths)   
    #Reshape, long to wide
    dists <- dcast(live,segmentID~obs,value.var='dist')
    names(dists)[2:ncol(dists)] <- paste0('dist.',names(dists)[2:ncol(dists)])

  #Add number of observations on each segment to survey data, including zeroes for segments with no observations
    
    surveys <- join(allSurveys,nDetects,by='segmentID',type='left')
    surveys$nDetects[is.na(surveys$nDetects)] <- 0
    
  #Add observation distances to survey data, then sort to bring surveys with observations to top
    
    surveys <- join(surveys,dists,by='segmentID',type='left')
    surveys <- surveys[with(surveys,order(-nDetects,recovID,season,segmentID)),] 

#-----------------------------------------------------------------------------------------------# 
# Save objects needed for analysis into an .Rdata file 
#-----------------------------------------------------------------------------------------------#   

  # save(g0obs,surveys,minYear,maxYear,Trunc,W,file='MDT_Data.Rdata')

    