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
    
    library(plyr)
    library(reshape2)
    library(raster)
    library(rgdal)

#-----------------------------------------------------------------------------------------------# 
# Import data
#-----------------------------------------------------------------------------------------------# 
  
  # Distance-sampling data
    
    # Survey data (e.g., segment names, lengths)
    allSurveys_orig <- read.csv('TortoiseData/Tort_Surveys.csv',
                                header=TRUE,na.strings=c("NA",""),strip.white=TRUE)
    allSurveys_1921 <- read.csv('TortoiseData/Tort_Surveys_2019-2021.csv',
                                header=TRUE,na.strings=c("NA",""),strip.white=TRUE)
    
    # Detection data for live tortoises
    live_orig <- read.csv('TortoiseData/Tort_Obs_Live.csv',header=TRUE,na.strings=c("NA",""),strip.white=TRUE)
    live_1921 <- read.csv('TortoiseData/Tort_Obs_Live_2019-2021.csv',header=TRUE,na.strings=c("NA",""),strip.white=TRUE)
    
    # Recovery units associated with established strata (or TCAs)
    ru <- read.csv('TortoiseData/RUs_TCAs.csv',header=TRUE,strip.white=TRUE)
    
  # Covariate data
    
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
    
    # Winter precipitation data (Oct-Mar preceding each survey, mm)
    ppt.surveys <- read.csv('Covariates/survey_mid_points_with_precip_DAYMET.csv',
                            header=TRUE,strip.white=TRUE)
    ppt.octmar <- raster('Covariates/30yr_oct_mar_prcp_mm.tif')          #Winter precip (Oct-Mar, mm), 30yr norms

  # Telemetry data (for estimating availability/g0)
    
    # Telemetry data, with binary "visible" variable (excluding RM site)
    g0obs_orig <- read.csv('TortoiseData/g0_Obs.csv',
                           header=TRUE,na.strings=c("NA",""),strip.white=TRUE) 
    g0obs_1921 <- read.csv('TortoiseData/g0_Obs_2019-2021.csv',
                           header=TRUE,na.strings=c("NA",""),strip.white=TRUE) 
    
    # Telemetry data from RM (teaching) site, 2012-2013 only
    g0obs_rm <- read.csv('TortoiseData/g0_Obs_RMonly.csv',
                         header=FALSE,na.strings=c("NA",""),strip.white=TRUE,
                         col.names=names(g0obs_orig))                                                    
    
    # Name, abbreviation, and location of each telemetry site
    g0sites <- read.csv('TortoiseData/g0_Sites.csv',
                        header=TRUE,na.strings=c("NA",""),strip.white=TRUE)
    
    # Winter precipitation data (Oct-Mar preceding each telemetry season, mm)
    ppt.g0 <-read.csv('Covariates/telemetry_sites_yearly_precip_DAYMET.csv',
                      header=TRUE,strip.white=T)
    
  # Shapefiles 
    
    # Recovery units
    rus <- readOGR(dsn='Covariates/Revised Recovery Units',layer='2011RecoveryUnits')
    rus <- spTransform(rus,crs(east))
    
    # TCAs
    tcas <- readOGR(dsn='Covariates/TCAs',layer='All_Strata')
    tcas <- spTransform(tcas,crs(east)) 

#-----------------------------------------------------------------------------------------------# 
# Establish start and end years
#-----------------------------------------------------------------------------------------------#  
  
    minYear <- 2001
    maxYear <- 2020    
    
#-----------------------------------------------------------------------------------------------# 
# Remove data from Upper Virgin River (UVR) recovery unit and a couple smaller regions
#-----------------------------------------------------------------------------------------------#  
  
  # Upper Virgin River
  
    #Remove data from survey data and live tortoise observation dataset
    allSurveys_orig <- allSurveys_orig[allSurveys_orig$stratum!='RC',]
    live_orig <- live_orig[live_orig$segmentID %in% allSurveys_orig$segmentID,]
    
    # Remove telemetry data from UVR (site = RC)
    g0obs_orig <- g0obs_orig[g0obs_orig$G0_site!='RC',]
    
  # Non_PI and Non_FK (few surveys in Searchlight and Red Mountain, respectively)

    allSurveys_orig <- allSurveys_orig[!allSurveys_orig$stratum %in% c('Non_PI','Non_FK'),]
    
#-----------------------------------------------------------------------------------------------# 
# Combine data through 2018 with 2019-2020 data
#-----------------------------------------------------------------------------------------------#  
  
  # Survey data
   
    allSurveys_orig <- allSurveys_orig[,c('year','stratum','transectID','segmentID','date_',
                                          'seg_length','group','start_lead','st_wp_sequ','end_wp_seq',
                                          'start_east','start_nort','end_eastin','end_northi')]
    allSurveys_1921 <- allSurveys_1921[,c('year_','stratum','transectID','segmentID','date_',
                                          'seg_length_m','group_','start_lead','st_wp_sequ','end_wp_sequ',
                                          'start_easting','start_northing','end_easting','end_northing')] 
    names(allSurveys_orig) <- names(allSurveys_1921) <- c('year','stratum','transectID','segmentID','date',
                                                          'seg_length','group','start_lead','st_wp_seq',
                                                          'end_wp_seq','start_easting','start_northing',
                                                          'end_easting','end_northing')
    allSurveys <- rbind(allSurveys_orig,allSurveys_1921)
    # Check:
    # summary(allSurveys[,c(1,5,6,11:14)])
    # summary(allSurveys_orig[,c(1,5,6,11:14)])

  # Live tortoise observations
    
    live_1921$mcl_greater_180 <- ifelse(live_1921$mcl_greater_180=='Yes','yes',
                                        ifelse(live_1921$mcl_greater_180=='No','no','unknown'))
    live <- rbind(live_orig,live_1921)

  # Telemetry data (also merging RM data with everything else)
    
    g0obs <- rbind(g0obs_orig,g0obs_rm,g0obs_1921)
    
  # For each dataset, only include data from selected years (2001-2020)
    
    allSurveys <- allSurveys[allSurveys$year %in% 2001:2020,]
    live <- live[live$year %in% 2001:2020,]
    g0obs <- g0obs[g0obs$year_ %in% 2001:2020,]

#-----------------------------------------------------------------------------------------------# 
# Calculate midpoint of each segment
#-----------------------------------------------------------------------------------------------#  
  
  allSurveys$mid_easting <- round((allSurveys$start_easting + allSurveys$end_easting)/2)
  allSurveys$mid_northing <- round((allSurveys$start_northing + allSurveys$end_northing)/2)
        
#-----------------------------------------------------------------------------------------------# 
# Add recovery unit to survey data, live observations
#-----------------------------------------------------------------------------------------------#     

  # Adding recovery unit to survey data
  
    allSurveys$recovID_incomplete <- ru$recovery_unit[match(allSurveys$stratum,ru$stratum)]
    # 314 NAs for strata that weren't in ru file (because they aren't associated with a single RU)
    
    count(allSurveys[,c('recovID_incomplete','stratum')])   
      # Colorado Desert: AG, CK, CM, FE, JT, PT, PV (each with with 1262 - 2318 segments)
      # Eastern Mojave: EV, IV (2175, 2176 segments)
      # Eastern Mojave: PN, PS (Pahrump North and South; 158,282 segments)
      # Northeastern Mojave: BD, CS, GB, MM (2282-3512 segments)
      # Northeastern Mojave: BD2, MM2, Non_CS (10-362 segments)
      # Western Mojave: FK, OR, SC (1993-3958 segments)
      # Western Mojave: EF, MC, Non_SC (Edwards AFB, Marine Corps AGCC 29 Palms; 8-190 segments)
      # No recovery unit (b/c in multiple): MN, MS, MP, NS, PT2 
        # (Lake Mean North and South, Mojave Nat Preserve, Newberry Springs; 15-102 segments)
    
    segmids <- allSurveys[,c('segmentID','recovID_incomplete','mid_easting','mid_northing')]
    seglocs_mid <- SpatialPointsDataFrame(coords=segmids[,3:4],proj4string=crs(elev),data=segmids)
    
    # Extract name of RU that segment mid point falls into (based on shapefile)
    segmids.new <- over(seglocs_mid,rus[,'Unit_Name'])  

    # Use recovery unit assignment based on polygon overlay where needed
    allSurveys <- cbind(allSurveys,segmids.new)
    allSurveys$recovID <- ifelse(is.na(allSurveys$recovID_incomplete),
                                 allSurveys$Unit_Name,allSurveys$recovID_incomplete)

    # Did a quick check and the new assignments all seem correct
    # Check that all segments in a transect have the same RU assignment
      # trans <- ddply(allSurveys,.(transectID),summarize,nSegs=length(unique(segmentID)),nRU=length(unique(recovID)))
      # summary(trans) #Yes, there's only one RU assignment for each transect
    
    allSurveys$recovID_incomplete <- NULL
    allSurveys$Unit_Name <- NULL

    live$recovID <- allSurveys$recovID[match(live$segmentID,allSurveys$segmentID)]
    
#-----------------------------------------------------------------------------------------------# 
# Format telemetry data
#-----------------------------------------------------------------------------------------------#  

  # Format date, year, and day-of-year
    
    g0obs$obsdate <- as.Date(g0obs$date_,format="%m/%d/%Y")  
    # Recreate year from date (and name it "season" to match precipitation data)
    g0obs$season  <- as.numeric(format(g0obs$obsdate,'%Y'))   
    g0obs$doy <- as.numeric(format(g0obs$obsdate,'%j'))

  # Create unique tortoise IDs because tort_num is not unique across sites
  
    g0obs$tortid <- paste(g0obs$G0_site,g0obs$tort_num,sep='_')  
    
  # Convert response (visible) to numeric values
  
    # count(g0obs$visible)
    g0obs <- g0obs[g0obs$visible!='Unk',] 
    g0obs$visible[g0obs$visible %in% c('N','No')]  <- 0
    g0obs$visible[g0obs$visible %in% c('Y','Yes')] <- 1
    g0obs$visible <- as.numeric(g0obs$visible)
  
  # Make sure that all sites in g0obs are in g0sites
  
    #Should return "character(0)" if all sites in g0obs are in g0sites
    unique(g0obs$G0_site[!g0obs$G0_site %in% g0sites$G0_site]) 
    
  # Add site UTMs, max temp normals, and Oct-Mar precipitation normals to g0obs

    g0sites <- g0sites[,c('G0_site','Easting','Northing')]
    g0locs <- SpatialPointsDataFrame(coords=g0sites[,2:3],proj4string=crs(ppt.octmar),data=g0sites)
    g0sites$temp.max <- extract(temp.max,g0locs)
    g0sites$ppt.norm <- extract(ppt.octmar,g0locs)
    g0obs <- join(g0obs,g0sites,by='G0_site',type='left')

  # Add yearly precipitation data to g0obs
    
    # Check that both sets of data include the same sites
    all.equal(sort(ppt.g0$g0site),sort(unique(g0obs$G0_site)))
    
    # Clean up yearly precipitation data and put in long form
    ppt.g0 <- ppt.g0[,c(1,5:ncol(ppt.g0))]
    names(ppt.g0) <- c('G0_site',paste('ppt',2001:2020,sep='.'))
    ppt.g0.long <- melt(ppt.g0,id.vars='G0_site',variable.name='ppt.yr',value.name='ppt')
    ppt.g0.long$season <- as.numeric(substr(ppt.g0.long$ppt.yr,5,8))
    
    # Add to g0obs
    g0obs <- join(g0obs,ppt.g0.long[,c('G0_site','ppt','season')],by=c('G0_site','season'),type='left')
    
  # Calculate annual % of 30-yr Oct-Mar precipitation norms 
    
    g0obs$pptOM.perc <- g0obs$ppt/g0obs$ppt.norm*100
    # summary(g0obs$pptOM.perc) #ranges from 7-315%

  # Keep only the columns we'll need and rename
    
    g0obs <- g0obs[,c('tortid','G0_site','Easting','Northing','temp.max','obsdate','season','doy','visible','ppt','ppt.norm','pptOM.perc')]
    names(g0obs)[c(2:4)] <- c('site','east','north')  

#-----------------------------------------------------------------------------------------------# 
# Format observations of live tortoises (each row in dataframe represents one tortoise)
#-----------------------------------------------------------------------------------------------# 

  # Round distances to 1 decimal place 
    
    live$perp_distance_m <- round(live$perp_distance_m,1) 

  # Fix indicator, so mcl_greater_180='yes' when mcl measurement >= 180 mm:
    
    live$mcl_greater_180[which(live$mcl_mm>=180)] <- 'yes'
    
  # Remove unnecessary variables (keeping mcl variables and sex so that we could subset data by either) and rename columns
      
    live <- live[,names(live) %in% c('segmentID','year_','perp_distance_m',
                                     'sex','mcl_mm','mcl_greater_180'),]
    names(live)[c(1,3)] <- c('season','dist')

  # Subset observations to exclude tortoises with MCL<180 
  # (or exclude if we don't know whether MCL>=180)

    live <- live[live$mcl_greater_180=='yes',]  
    
  # Truncate observations by detection distance
    
    # If truncating, set Trunc = W = truncation distance (Eguchi and Gerrodette 2009)
    # If not truncating, set W = maxDist (max. perpendicular distance). 
    #Note: need to change the likelihood in STAN model if we want to do this!
    
    Trunc <- W <- 20   
    # nrow(live) #5647
    live <- live[live$dist<=Trunc,]
    # nrow(live) #5392 (truncating 4.5% of observations)
    
    liveyr <- ddply(live,.(season),summarize,ndetects=length(dist))
    summary(liveyr) #69-542 torts detected (within trunction distance) each year
    
    live <- live[,c('segmentID','dist')]
    
#-----------------------------------------------------------------------------------------------# 
# Format survey data (each row in dataframe represents one segment)
#-----------------------------------------------------------------------------------------------#  
 
  # Remove unnecessary columns
    
    allSurveys <- allSurveys[,!names(allSurveys) %in% c('st_wp_seq','end_wp_seq')]
    
  # Rename year column 
    
    names(allSurveys)[names(allSurveys)=='year'] <- 'season'

  # Format date, calculate day-of-year
    
    allSurveys$sdate <- as.Date(allSurveys$date,format='%m/%d/%Y')
    allSurveys$doy <- as.numeric(format(allSurveys$sdate,'%j'))
    allSurveys$date <- NULL

#-----------------------------------------------------------------------------------------------# 
# Add Oct-Mar precipitation data (ppt.surveys) to the survey data
#-----------------------------------------------------------------------------------------------#    

  # First, check that all segments in allSurveys are in the ppt.surveys dataframe
    
    # Should return "character(0)" if all segments in allSurveys are in ppt.surveys
    unique(allSurveys$segmentID[!allSurveys$segmentID %in% ppt.surveys$segmentID])

  # Add annual precipitation data to survey data
    
    ppt.surveys <- ppt.surveys[,c('prev_6mo_precip','segmentID')]
    names(ppt.surveys)[1] <- 'ppt'
    allSurveys$ppt <- ppt.surveys$ppt[match(allSurveys$segmentID,ppt.surveys$segmentID)]
    
  # Add Oct-Mar precipitation normals to survey data
    
    mids <- allSurveys[,c('segmentID','mid_easting','mid_northing')]
    midlocs <- SpatialPointsDataFrame(coords=mids[,2:3],proj4string=crs(elev),data=mids)
    allSurveys$ppt.norm <- extract(ppt.octmar,midlocs)
    
  # Calculate annual % of 30-yr Oct-Mar precipitation norms 
    
    allSurveys$pptOM.perc <- allSurveys$ppt/allSurveys$ppt.norm*100

#-----------------------------------------------------------------------------------------------# 
# Attach spatial covariate data to allSurveys
#-----------------------------------------------------------------------------------------------#   
    
  # The wash layer is missing some values in the north-central part of the study area
  # Will impute the mean value
  
    wash.impute <- wash
    wash.impute[is.na(wash.impute)] <- mean(values(wash),na.rm=TRUE)
    
  # Extract values for each segment based on start location
    
    mids$precip.s <- extract(precip.s,midlocs)
    mids$precip.w <- extract(precip.w,midlocs)
    mids$temp.max <- extract(temp.max,midlocs)
    mids$elev <- extract(elev,midlocs)
    mids$slope <- extract(slope,midlocs)
    mids$aspect.n <- extract(north,midlocs)
    mids$aspect.e <- extract(east,midlocs)
    mids$rough <- extract(rough,midlocs)
    mids$wash <- extract(wash.impute,midlocs)
    mids$bedrock <- extract(bedrock,midlocs)
    mids$soilbulk <- extract(soilbulk,midlocs)
    mids$prock <- extract(prock,midlocs)
    mids$veg.p <- extract(veg.p,midlocs)
    mids$veg.a <- extract(veg.a,midlocs)
    mids$road.m <- extract(road.m,midlocs)
    mids$road.a <- extract(road.a,midlocs)
    
    # Check that there are no missing values:
    summary(mids)

  # Merge covariate data with survey data

    allSurveys <- join(allSurveys,mids[,c(1,4:ncol(mids))],by='segmentID',type='left')
    
#-----------------------------------------------------------------------------------------------# 
# Merge survey and observation data
#-----------------------------------------------------------------------------------------------# 

  # Count the number of observations per segment
    
    nDetects <- as.data.frame(table(live$segmentID)) 
    names(nDetects) <- c('segmentID','nDetects')
    
  # Reshape detection distances from long to wide (all observations for a segment on one row)
    
    live <- live[with(live,order(segmentID,dist)),]      
    #Add observation number (1:nDetects) for each detection on each segment
    live$obs <- sequence(rle(as.character(live$segmentID))$lengths)   
    #Reshape, long to wide
    dists <- dcast(live,segmentID~obs,value.var='dist')
    names(dists)[2:ncol(dists)] <- paste0('dist.',names(dists)[2:ncol(dists)])

  # Add number of observations on each segment to survey data, 
  # including zeroes for segments with no observations
    
    surveys <- join(allSurveys,nDetects,by='segmentID',type='left')
    surveys$nDetects[is.na(surveys$nDetects)] <- 0
    
  # Add observation distances to survey data, then sort to bring surveys with observations to top
    
    surveys <- join(surveys,dists,by='segmentID',type='left')
    surveys <- surveys[with(surveys,order(-nDetects,recovID,season,segmentID)),] 

#-----------------------------------------------------------------------------------------------# 
# Clean up stratum assignments in survey data
#-----------------------------------------------------------------------------------------------# 

  count(surveys$stratum)
  
  # Stratum is not used as a covariate in any part of the model, but it could be useful to assess survey effort
  # Create new assignment so I keep the original
  
    surveys$stratum_new <- surveys$stratum
    surveys$stratum_new[surveys$stratum_new=='BD2']    <- 'BD'
    surveys$stratum_new[surveys$stratum_new=='MM2']    <- 'MM'
    surveys$stratum_new[surveys$stratum_new=='Non_CS'] <- 'CS'
    surveys$stratum_new[surveys$stratum_new=='Non_SC'] <- 'SC'
    surveys$stratum_new[surveys$stratum_new=='PT2']    <- 'PT'
    
    # Consider EF (Edwards AFB) part of Fremont-Kramer
    surveys$stratum_new[surveys$stratum_new=='EF'] <- 'FK'    
    
    # Consider MP (Mojave Preserve outside of IV) part of IV
    surveys$stratum_new[surveys$stratum_new=='MP'] <- 'IV'    
    
    # Combine Pahrump North and South (only surveyed in 2008)
    surveys$stratum_new[surveys$stratum_new %in% c('PN','PS')] <- 'PA'

    # Combine Lake Mead North and South (NPS; only surveys in 2001, 2005)
    surveys$stratum_new[surveys$stratum_new %in% c('MN','MS')] <- 'ME'

    # Keep MC (MCAGCC, 29 Palms) on it's own
    # Keep NS (Newberry Springs [designated wilderness]) on it's own
    
    # Now we're left with 20 "strata"
    # sty <- ddply(surveys,.(stratum_new,recovID,season),summarize,ntransects=length(unique(transectID)))
    # sty.w <- dcast(sty,stratum_new + recovID ~ season)

#-----------------------------------------------------------------------------------------------# 
# Summarizing survey dates
#-----------------------------------------------------------------------------------------------#       
  
  surveys$mon <- as.numeric(format(surveys$sdate,'%m'))
  surveys$day <- as.numeric(format(surveys$sdate,'%d'))

  dates <- count(surveys[,c('doy','mon','day')]) #65 (5-Mar) to 179 (28-Jun)
    head(dates)
    tail(dates)

  quantile(surveys$doy,c(0.05,0.10,0.25,0.5,0.75,0.9,0.95))
    # 90% of surveys done with doy = 79 (20-Mar) - 147 (27-May)

#-----------------------------------------------------------------------------------------------# 
# Summarizing survey effort
#-----------------------------------------------------------------------------------------------# 
  
  # Number of transects surveyed each year
    seasons <- ddply(surveys,.(season),summarize,nru=length(unique(recovID)),
                     nstratum=length(unique(stratum_new)),ntran=length(unique(transectID)),
                     nseg=length(unique(segmentID)),nkm=sum(seg_length)/1000)
    seasons[seasons$season<2004,]
    # 2001-2003 (1728, 1068, 990)
    summary(seasons[seasons$season>2003,])
    # 2004-2020 (173-894)
    
  # TCAs/strata:
    length(unique(surveys$stratum_new)) #20 (including 16 in shapefile + PA + ME + MC + NS)
    stratayr <- ddply(surveys,.(season,stratum_new),summarize,ntran=length(unique(transectID)),
                      nseg=length(unique(segmentID)))
    summary(stratayr)
    summary(stratayr[stratayr$season %in% 2001:2003,]) # 12-299 transects per strata
    summary(stratayr[stratayr$season %in% 2004:2020,]) # 5-198 transects per strata

  # Number of detections per segment
    summary(surveys$nDetects) # mean = 0.1394
    count(surveys$nDetects)

  # Number of detections per km
    sum(surveys$nDetects)/sum(surveys$seg_length/1000) # 0.0482

  # Segment lengths 
    ll <- ddply(surveys,.(season),summarize,min=min(seg_length),max=max(seg_length),
                mean=round(mean(seg_length)),med=median(seg_length))  
    # Median length in 2001: 1600 m
    # Median lengths in 2002-2003: 3996 m
    # Median length in 2004: 2505 m
    # Median lengths in 2005-2020: 2995-3003 m
    
#-----------------------------------------------------------------------------------------------# 
# Summarizing radiotelemetry data
#-----------------------------------------------------------------------------------------------#  

  g0seas <- ddply(g0obs,.(season),summarize,nsite=length(unique(site)),ntort=length(unique(tortid)))
  g0seassite <- ddply(g0obs,.(season,site),summarize,ntort=length(unique(tortid)))
  summary(g0seassite) # 4-23 tortoises per site/yr

  # Number of radiomarked tortoises
    length(unique(g0obs$tortid)) # 516

  # Mean availability
    mean(g0obs$visible) # 0.7722
    sd(g0obs$visible)/sqrt(nrow(g0obs)) # 0.00175
    
#-----------------------------------------------------------------------------------------------# 
# Save objects needed for analysis into an .Rdata file 
#-----------------------------------------------------------------------------------------------#   

  # save(g0obs,surveys,minYear,maxYear,Trunc,W,file='MDT_Data.Rdata'
