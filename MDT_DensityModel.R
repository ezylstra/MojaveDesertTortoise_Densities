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
    library(rstan)
    rstan_options(auto_write = TRUE)
    rstan_options(javascript = FALSE)
    options(mc.cores = parallel::detectCores())

#-----------------------------------------------------------------------------------------------# 
# Import data
#-----------------------------------------------------------------------------------------------# 
  
  load('MDT_Data.Rdata') 

#-----------------------------------------------------------------------------------------------# 
# Format covariates associated with telemetry observations (explaining variation in availability)
#-----------------------------------------------------------------------------------------------#   

  # Standardize potential covariates (z-scores) and save as separate objects
    
    telem.easting.mn  <- mean(g0obs$east)
    telem.northing.mn <- mean(g0obs$north)
    telem.doy.mn      <- mean(g0obs$doy)
    telem.precip.mn   <- mean(g0obs$pptOM.perc)
    telem.easting.sd  <- sd(g0obs$east)
    telem.northing.sd <- sd(g0obs$north)
    telem.doy.sd      <- sd(g0obs$doy)
    telem.precip.sd   <- sd(g0obs$pptOM.perc)
    
    telem.easting  <- (g0obs$east - telem.easting.mn)/telem.easting.sd
    telem.northing <- (g0obs$north - telem.northing.mn)/telem.northing.sd
    telem.doy      <- (g0obs$doy - telem.doy.mn)/telem.doy.sd
    telem.precip   <- (g0obs$pptOM.perc - telem.precip.mn)/telem.precip.sd
    
  # Create quadratics
    
    telem.easting2   <- telem.easting * telem.easting
    telem.northing2  <- telem.northing * telem.northing
    telem.doy2       <- telem.doy * telem.doy
    telem.precip2    <- telem.precip * telem.precip
  
#-----------------------------------------------------------------------------------------------# 
# Create indices for year, recovery units, transects
#-----------------------------------------------------------------------------------------------#  

    year <- as.integer(surveys$season-minYear + 1)   # Set minYear as year one
    surveys$recovID.ind <- as.integer(as.factor(surveys$recovID)) # CD = 1; EM = 2; NEM = 3; WM = 4  
    recovID <- surveys$recovID.ind   
    surveys$transectID.ind <- as.integer(as.factor(surveys$transectID))
    transectID <- surveys$transectID.ind
    
#-----------------------------------------------------------------------------------------------# 
# Format covariates associated with distance-sampling SEGMENTS
#-----------------------------------------------------------------------------------------------#
    
  # Transform and standardize potential covariates, then save as separate objects
    # Only transforming variables when necessary to keep scaled values in a reasonable range.
    # Log transforming assumes a change in lower values has a bigger effect than a change at higher values

    precip.s <- scale(surveys$precip.s)        #Summer precip (May-Oct; dry season), 30-year normals (mm)
    precip.w <- scale(surveys$precip.w)        #Winter precip (Nov-Apr; wet season), 30-year normals (mm)
    temp.max <- scale(surveys$temp.max)        #Temperature in warmest month, 30-year normals (degC)
    elev <- scale(surveys$elev)                #Elevation (m)
    slope <- scale(log(surveys$slope))         #Slope (degree)
    aspect.n <- scale(surveys$aspect.n)        #Northness (1 = north; -1 = south)
    aspect.e <- scale(surveys$aspect.e)        #Eastness (1 = east; -1 = west)
    rough <- scale(log(surveys$rough))         #Average surface roughness (ratio of area to planimetric area)
    wash <- scale(log(surveys$wash))           #Proportion of grid cell classified as "wash"
    bedrock <- scale(log(surveys$bedrock))     #Depth to bedrock (cm)
    soilbulk <- scale(surveys$soilbulk)        #Soil bulk density (kg/m3)
    prock <- scale(surveys$prock)              #Percent rock >254 mm
    veg.p <- scale(log(surveys$veg.p))         #Perennial vegetation (NDVI value in very dry year [2019])
    veg.a <- scale(surveys$veg.a)              #Annual plant (forage) potential (comparison of MODIS values in dry, wet years)
    road.m <- scale(log(surveys$road.m))       #Distance to minor road (m)
    road.a <- scale(log(surveys$road.a))       #Distance to any road (m)
    annprecip <- scale(surveys$pptOM.perc)     #Precipitation during previous winter (Oct-Mar, % 30yr norms) -- time varying
    
    lateyrs <- 1*(surveys$season>2003)         #Indicator for all years after 2003 

    # Create quadratic for maximum temperature
    temp.max2 <- temp.max * temp.max

  # Evaluate potential correlations (highlighting those with abs(r) >= 0.5):
    
    round(cor(cbind(precip.w,precip.s,temp.max,
                    slope,rough,aspect.n,aspect.e,
                    wash,bedrock,
                    veg.p,veg.a,
                    road.a)),3)

    round(cor(cbind(precip.w,precip.s,temp.max,rough,aspect.n,aspect.e,
                    wash,bedrock,veg.a,road.a)),3)  
    # With this combination, most extreme r = -0.498 (rough-bedrock)
    
  # Covariates that can potentially be used to predict availability on a given survey
    # Note: values drawn from surveys dataframe, but standardized based on mean and SD from telemetry data
    
    g0.doy <- (surveys$doy - telem.doy.mn)/telem.doy.sd
    g0.easting <- (surveys$mid_easting - telem.easting.mn)/telem.easting.sd  ### this used to be start_e
    g0.northing <- (surveys$mid_northing - telem.northing.mn)/telem.northing.sd  ### this used to be start_n
    g0.precip <- (surveys$pptOM.perc - telem.precip.mn)/telem.precip.sd
    
    g0.doy2 <- g0.doy * g0.doy
    g0.easting2 <- g0.easting * g0.easting
    g0.northing2 <- g0.northing * g0.northing
    g0.precip2 <- g0.precip * g0.precip

#-----------------------------------------------------------------------------------------------# 
# Select covariates for each "sub-model"
#-----------------------------------------------------------------------------------------------#
    
  # Choose abundance covariates at the segment level [time-invariant]
  # cov_lam[n_segments, n_cov_lam] 
    
    cov_lam <- data.frame(precip.w, precip.s, temp.max, temp.max2, rough, aspect.e, 
                          aspect.n, wash, bedrock, veg.a, road.a)
    cov_lam <- cov_lam[,sort(names(cov_lam)),drop=F]
    cov_lam <- as.matrix(cov_lam)
    n_cov_lam <- ncol(cov_lam)
    
  # Choose (segment-level) covariates for sigma, the SD in a half-normal detection function  
  # cov_detect[n_segments, n_cov_detect] 
    
    cov_detect <- data.frame(rough, annprecip, veg.p, lateyrs)
    cov_detect <- cov_detect[,sort(names(cov_detect)),drop=F]
    cov_detect <- as.matrix(cov_detect)
    n_cov_detect <- ncol(cov_detect)

  # Choose covariates affecting above-ground activity of radiomarked tortoises
  # cov_telem[n_telem, n_cov_telem]
    # Note: the covariates here should match those in cov_g0 (below)
    
    cov_telem <- data.frame(telem.easting, telem.easting2, telem.northing, telem.northing2, 
                            telem.doy, telem.doy2, telem.precip, telem.precip2)
    cov_telem <- cov_telem[,sort(names(cov_telem)),drop=F]
    cov_telem <- as.matrix(cov_telem)
    n_cov_telem <- ncol(cov_telem)

  # Choose covariates affecting above-ground activity of tortoises on distance sampling segments
  # cov_g0[n_segments, n_cov_telem]
    # Note: the covariates here should match those in cov_telem (above)

    cov_g0 <- data.frame(g0.easting, g0.easting2, g0.northing, g0.northing2,
                         g0.doy, g0.doy2, g0.precip, g0.precip2)
    cov_g0 <- cov_g0[,sort(names(cov_g0)),drop=F]
    cov_g0 <- as.matrix(cov_g0)
    
  # Calculate means of scaled, abundance-related covariates within each recovery unit
    
    segsc <- data.frame(recovID=as.numeric(as.factor(surveys$recovID)),
                        precip.w,precip.s,temp.max,elev,slope,rough,aspect.n,aspect.e,  
                        wash,bedrock,prock,soilbulk,veg.p,veg.a,road.m,road.a)
    segCov_RUmeans <- ddply(segsc,.(recovID),summarize, 
                            precip.w=mean(precip.w),precip.s=mean(precip.s),
                            temp.max=mean(temp.max), elev=mean(elev),slope=mean(slope),
                            rough=mean(rough), aspect.n=mean(aspect.n), aspect.e=mean(aspect.e),  
                            wash=mean(wash),bedrock=mean(bedrock),prock=mean(prock),
                            soilbulk=mean(soilbulk),veg.p=mean(veg.p),veg.a=mean(veg.a),
                            road.m=mean(road.m),road.a=mean(road.a))
    # Add quadratic for temp.max2 
    segCov_RUmeans$temp.max2 <- segCov_RUmeans$temp.max * segCov_RUmeans$temp.max
    
  # Create dataframe with just those covariates selected for cov_lam above
    
    RUmeans <- as.matrix(segCov_RUmeans[,names(segCov_RUmeans) %in% colnames(cov_lam)])
    RUmeans <- RUmeans[,sort(colnames(cov_lam))]     
    
#-----------------------------------------------------------------------------------------------# 
# Organize and summarize data for STAN
#-----------------------------------------------------------------------------------------------# 

  # Telemetry data
    
    y_telem <- g0obs$visible                          #Binary response variable (visible: 1/0)
    rmtortID  <- as.integer(as.factor(g0obs$tortid))  #Vector with index that uniquely identifies each tortoise
    n_telem_obs <- length(y_telem)                    #Number of telemetry observations
    n_rmtorts <- max(rmtortID)                        #Number of unique tortoises tracked

   # Distance-sampling data
    
    n_years <- maxYear - minYear + 1                 #Number of years included in trend/study 
    n_recov <- max(recovID)                          #Number of recovery units
    n_transects <- max(transectID)                   #Number of transects
    n_segments <- nrow(surveys)                      #Number of segments
    n_segments1 <- sum(surveys$nDetects>0)           #Number of segments with at least one detection
    yeartrend <- 1:n_years                           #Sequence of years, with first year = 1
    
    L <- surveys$seg_length                          #Vector with length of each segment surveyed (in meters)
    n <- surveys$nDetects                            #Vector with number of tortoise detections per segment
    
    # Creating matrix of distances for only those segments with >=1 detection
    # Converting NAs to 0s for STAN (won't matter since we'll loop over the non-zero values in each row)
    y <- as.matrix(surveys[1:n_segments1,grep('dist',names(surveys))])    
    y[is.na(y)] <- 0                                                    
  
  # Create dataframe that lists every year-recovery unit combination when surveys were done
    
    ru.yr <- unique(surveys[,c('recovID.ind','season')])
    ru.yr$yr.ind <- as.integer(ru.yr$season - minYear + 1) 
    pred.df <- data.frame(ru=ru.yr$recovID.ind,yr=ru.yr$yr.ind)
    pred.df <- pred.df[with(pred.df,order(ru,yr)),]
    n_preds <- nrow(pred.df)  #Number of year-RU combinations when surveys were done

#-----------------------------------------------------------------------------------------------# 
# Bundle everything for STAN
#-----------------------------------------------------------------------------------------------# 

  # MCMC parameters

    ni <- 3000    #No. iterations (including warmup) 
    nb <- 500     #No. burn-in iterations to discard (ie, warmup)
    nt <- 1       #Thin rate
    nc <- 3       #No. chains
    
  # Data
    
    tort_data <- list(
      n_recov=n_recov,
      n_transects=n_transects, 
      n_segments=n_segments, 
      n_segments1=n_segments1,
      n_years=n_years,
      n_telem_obs=n_telem_obs,
      n_rmtorts=n_rmtorts,
      transectID=transectID,
      recovID=recovID,
      year=year,
      yeartrend=yeartrend,
      lateyrs=lateyrs,
      rmtortID=rmtortID,
      cov_lam=cov_lam,
      cov_detect=cov_detect,
      cov_g0=cov_g0,
      cov_telem=cov_telem,
      n_cov_lam=n_cov_lam,
      n_cov_detect=n_cov_detect,
      n_cov_telem=n_cov_telem,
      y_telem=y_telem,
      n=n, 
      y=y,
      L=L, 
      W=W,
      T=Trunc,
      RUmeans=RUmeans,
      n_preds=n_preds,
      recov_pred=as.vector(pred.df$ru),
      yr_pred=as.vector(pred.df$yr))
      
  # Initial values
    
    set.seed(1)
    inits <- lapply(1:nc, function (i)
             list(beta_lam=runif(n_cov_lam,-0.5,0.5),
                  beta_detect=runif(n_cov_detect,-0.5,0.5),
                  beta_g0=runif(n_cov_telem,-0.5,0.5),
                  mu_recov=runif(n_recov,1,5),
                  trend=runif(n_recov,-1,1),
                  mu_sigma=runif(1,1,10),
                  mu_g0=runif(1,0.7,1),
                  sd_recov=runif(1,0,1),
                  sd_tran=runif(1,0,1),
                  sd_seg=runif(1,0,1),
                  sd_detect=runif(1,0,1),
                  sd_g0=runif(1,0,1),
                  recov_l_raw=matrix(runif(n_recov*n_years,-1,1),nrow=n_recov,ncol=n_years),
                  tran_l_raw=runif(n_transects,-1,1),
                  seg_l_raw=runif(n_segments,-1,1),
                  detect_raw=runif(n_segments,-1,1),
                  g0_raw=runif(n_rmtorts,-1,1)))

  # Parameters to monitor
    
    params <- c('beta_lam','beta_detect','beta_g0',
                'mu_recov','trend','mu_sigma','mu_g0',
                'sd_recov','sd_tran','sd_seg','sd_detect','sd_g0',
                'predTrend','D_survyrs','fit_obs','fit_new','recov_l_raw')
  
#-----------------------------------------------------------------------------------------------# 
# Call STAN from R
#-----------------------------------------------------------------------------------------------# 
  
  start.time <- Sys.time()  
  out <- stan('MDT_DensityModel.stan',
              control=list(adapt_delta=0.99),
              data=tort_data, init=inits, pars=params,
              chains=nc, iter=ni, warmup=nb, thin=nt,
              seed=1, cores=3, open_progress=FALSE)
  end.time <- Sys.time()
  
  # Save objects for other post-processing
  
    save(out,pred.df,surveys,start.time,end.time,
         cov_lam,cov_detect,cov_g0,cov_telem,
         file="STANFit_MDT_Density_2001-2020.Rdata")
  
  # Summary of model run  
    
    print(out,digits=2) 
    
#-----------------------------------------------------------------------------------------------# 
# Simple summaries of posterior distributions (with associated names of covariates)
#-----------------------------------------------------------------------------------------------# 

  posterior <- as.matrix(out)
  n_cov_lam <- ncol(cov_lam)
  n_cov_detect <- ncol(cov_detect)
  n_cov_g0 <- ncol(cov_g0)

  # Estimates of covariates (abundance, detection, and availability models)
    ncovar <- n_cov_lam + n_cov_detect + n_cov_g0
    stan.name <- rownames(summary(out)$summary)[1:ncovar]
    submodel <- c(rep('Lambda',n_cov_lam),rep('Detection',n_cov_detect),rep('Availability',n_cov_g0))
    covariate <- c(colnames(cov_lam),colnames(cov_detect),colnames(cov_g0))
    meansd <- round(as.data.frame(summary(out)$summary[1:ncovar,c(1,3)],row.names=F),2)
    quants <- round(as.data.frame(summary(out)$summary[1:ncovar,c(4,6,8)],row.names=F),2)
    names(quants) <- c('LCL','Median','UCL')
    cov.df <- cbind(stan.name,submodel,covariate,meansd,quants)
    cov.df$exclude0 <- ifelse(cov.df$LCL>0 & cov.df$UCL>0 | cov.df$LCL<0 & cov.df$UCL<0,'yes','no')
    cov.df

  # Estimates of log-linear trends
    trendindex <- grep('trend',rownames(summary(out)$summary))
    stan.name <- rownames(summary(out)$summary)[trendindex]
    trend <- paste('trend',c('CD','EM','NM','WM'),sep='-')
    trend.meansd <- round(as.data.frame(summary(out)$summary[trendindex,c(1,3)],row.names=F),3)
    trend.quants <- round(as.data.frame(summary(out)$summary[trendindex,c(4,6,8)],row.names=F),3)
    names(trend.quants) <- c('LCL','Median','UCL')
    trend.df <- cbind(stan.name,trend,trend.meansd,trend.quants)
    trend.df$exclude0 <- ifelse(trend.df$LCL>0 & trend.df$UCL>0 | trend.df$LCL<0 & trend.df$UCL<0,'yes','no')
    trend.df  

#-----------------------------------------------------------------------------------------------# 
# Assessing model fit
#-----------------------------------------------------------------------------------------------#  
    
    min.fit <- min(c(posterior[,'fit_obs'],posterior[,'fit_new']))-100
    max.fit <- max(c(posterior[,'fit_obs'],posterior[,'fit_new']))+100

    par(mfrow=c(1,1),mar=c(4,4,1,1))
    plot(posterior[,'fit_obs'],posterior[,'fit_new'],pch=19,cex=0.5,xlim=c(min.fit,max.fit),ylim=c(min.fit,max.fit))
    abline(a=0,b=1)
    (p.fit <- mean(posterior[,'fit_new'] > posterior[,'fit_obs']))  

#-----------------------------------------------------------------------------------------------# 
# Plot trends in each recovery unit
#-----------------------------------------------------------------------------------------------# 
  
  # Extract and format density estimates for RUs in each year they were surveyed
    
    densities <- posterior[,grep('D_',colnames(posterior))]
    pred.df$D.mn <- apply(densities,2,mean)
    pred.df$D.md <- apply(densities,2,median)
    pred.df$D.l <- apply(densities,2,quantile,probs=0.025)
    pred.df$D.u <- apply(densities,2,quantile,probs=0.975)
    pred.df$year <- pred.df$yr+2000
    dens <- pred.df    

  # Extract and format trend-based density estimates (every year and RU)
    
    pred.trend <- posterior[,grep('predTrend',colnames(posterior))]
    n_recov <- max(pred.df$ru)
    n_year <- max(pred.df$yr)
    trends <- data.frame(ru=rep(1:n_recov,n_year),year=rep(2001:(2001+n_year-1),each=n_recov),
                         t.mn=apply(pred.trend,2,mean),t.md=apply(pred.trend,2,median),
                         t.l=apply(pred.trend,2,quantile,probs=0.025),
                         t.u=apply(pred.trend,2,quantile,probs=0.975))

  # Create figure

    cexfix <- 0.8
    
    # pdf('Figures/Density_Trends.pdf', width = 3.15, height = 4.35, pointsize = 7)
    par(mfrow=c(n_recov,1),mar=c(0.5,1.0,0.5,1.0)+0.1,oma=c(2.5,2.5,0,0))
    # Colorado Desert
    plot(t.md~year,data=trends[trends$ru==1,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,10),bty='n',xaxs="i",yaxs='i',lwd = 0.5)
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F, lwd = 0.5)
      axis(1,at=seq(2002,2018,by=4),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0),lwd = 0.5,cex.axis = 1.2)
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F,lwd = 0.5)
      axis(2,at=seq(0,9,by=3),labels=seq(0,9,by=3),tcl=-0.25,las=1,mgp=c(1.5,0.5,0),lwd = 0.5, 
           cex.axis = 1.2)     
      lines(t.l~year,data=trends[trends$ru==1,],lty=2,lwd=0.5)
      lines(t.u~year,data=trends[trends$ru==1,],lty=2,lwd=0.5)
      arrows(x0=dens$year[dens$ru==1],y0=dens$D.l[dens$ru==1],x1=dens$year[dens$ru==1],y1=dens$D.u[dens$ru==1],
           col='black',length=0, lwd = 0.5)
      points(D.md~year,data=dens[dens$ru==1,],pch=19,cex=0.9)
      # mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      text(x=2020.2,y=9,adj=c(1,0),'Colorado Desert',cex=1.2)
    # Eastern Mojave
    plot(t.md~year,data=trends[trends$ru==2,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,10),bty='n',xaxs="i",yaxs='i',lwd = 0.5)
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F, lwd = 0.5)
      axis(1,at=seq(2002,2018,by=4),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0),lwd = 0.5,cex.axis = 1.2)
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F,lwd = 0.5)
      axis(2,at=seq(0,9,by=3),labels=seq(0,9,by=3),tcl=-0.25,las=1,mgp=c(1.5,0.5,0),lwd = 0.5, 
           cex.axis = 1.2)       
      lines(t.l~year,data=trends[trends$ru==2,],lty=2,lwd=0.5)
      lines(t.u~year,data=trends[trends$ru==2,],lty=2,lwd=0.5)
      arrows(x0=dens$year[dens$ru==2],y0=dens$D.l[dens$ru==2],x1=dens$year[dens$ru==2],y1=dens$D.u[dens$ru==2],
           col='darkgray',length=0,lwd=0.5)
      points(D.md~year,data=dens[dens$ru==2,],pch=19,cex=0.9)
      #mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      text(x=2020.2,y=9,adj=c(1,0),'Eastern Mojave',cex=1.2)
    # Northeastern Mojave
    plot(t.md~year,data=trends[trends$ru==3,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,10),bty='n',xaxs="i",yaxs='i',lwd = 0.5)
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F, lwd = 0.5)
      axis(1,at=seq(2002,2018,by=4),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0),lwd = 0.5,cex.axis = 1.2)
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F,lwd = 0.5)
      axis(2,at=seq(0,9,by=3),labels=seq(0,9,by=3),tcl=-0.25,las=1,mgp=c(1.5,0.5,0),lwd = 0.5, 
           cex.axis = 1.2)     
      lines(t.l~year,data=trends[trends$ru==3,],lty=2,lwd=0.5)
      lines(t.u~year,data=trends[trends$ru==3,],lty=2,lwd=0.5)
      arrows(x0=dens$year[dens$ru==3],y0=dens$D.l[dens$ru==3],x1=dens$year[dens$ru==3],y1=dens$D.u[dens$ru==3],
           col='darkgray',length=0,lwd=0.5)
      points(D.md~year,data=dens[dens$ru==3,],pch=19,cex=0.9)
      #mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      text(x=2020.2,y=9,adj=c(1,0),'Northeastern Mojave',cex=1.2)
    # Western Mojave
    plot(t.md~year,data=trends[trends$ru==4,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,10),bty='n',xaxs="i",yaxs='i',lwd = 0.5)
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F, lwd = 0.5)
      axis(1,at=seq(2002,2018,by=4),labels=seq(2002,2018,by=4),tcl=-0.25,mgp=c(1.5,0.5,0),lwd = 0.5,cex.axis = 1.2)
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F,lwd = 0.5)
      axis(2,at=seq(0,9,by=3),labels=seq(0,9,by=3),tcl=-0.25,las=1,mgp=c(1.5,0.5,0),lwd = 0.5, 
           cex.axis = 1.2)    
      lines(t.l~year,data=trends[trends$ru==4,],lty=2,lwd=0.5)
      lines(t.u~year,data=trends[trends$ru==4,],lty=2,lwd=0.5)
      arrows(x0=dens$year[dens$ru==4],y0=dens$D.l[dens$ru==4],x1=dens$year[dens$ru==4],y1=dens$D.u[dens$ru==4],
           col='darkgray',length=0,lwd=0.5)
      points(D.md~year,data=dens[dens$ru==4,],pch=19,cex=0.9)
      mtext('Density (tortoises/sq.km)',side=2,las=0,line=1.2,outer = TRUE)
      mtext('Year',side=1,line=2)
      text(x=2020.2,y=9,adj=c(1,0),'Western Mojave',cex=1.2)
    # dev.off()

#-----------------------------------------------------------------------------------------------# 
# Calculating marginal effects of covariates on predicted availability 
#-----------------------------------------------------------------------------------------------# 

  easting.mn <- 653821.2
  easting.sd <- 70420.52
  northing.mn <- 3923048
  northing.sd <- 112399.2
  doy.mn <- 113.854
  doy.sd <- 18.93502
  ppt.mn <- 99.01025
  ppt.sd <- 61.5989
  
  # Value of covariate for which predicted avilability is at a maximum:
    # maxima for y = ax^2 + bx + c: -b/2a
	
  	max.easting.z <- -cov.df$mean[cov.df$covariate=='g0.easting']/(2*cov.df$mean[cov.df$covariate=='g0.easting2'])
  	max.easting <- max.easting.z * easting.sd + easting.mn 
  	
  	max.doy.z <- -cov.df$mean[cov.df$covariate=='g0.doy']/(2*cov.df$mean[cov.df$covariate=='g0.doy2'])
  	max.doy <- max.doy.z * doy.sd + doy.mn 
  	
  	max.ppt.z <- -cov.df$mean[cov.df$covariate=='g0.precip']/(2*cov.df$mean[cov.df$covariate=='g0.precip2'])
    max.ppt <- max.ppt.z * ppt.sd + ppt.mn
   
  # Predicted availability with rainfall = max.ppt or equal to normal (at peak date):
    
    g0site.locs.z <- unique(cov_telem[,c('telem.easting','telem.easting2','telem.northing','telem.northing2')])
    beta_g0 <- posterior[,grep('beta_g0',colnames(posterior))]
    logita0 <- log(posterior[,'mu_g0']/(1-posterior[,'mu_g0']))
    
    # Mean predicted value on 19-April at each of the 22 sites with maximum winter precipitation
    X <- data.frame(doy=max.doy.z,doy2=max.doy.z*max.doy.z,
                    east=g0site.locs.z[,'telem.easting'],east2=g0site.locs.z[,'telem.easting2'],
                    north=g0site.locs.z[,'telem.northing'],north2=g0site.locs.z[,'telem.northing2'],
                    ppt=max.ppt.z,ppt2=max.ppt.z*max.ppt.z)
    logita <- matrix(1,nrow=22,ncol=1) %*% t(logita0) + as.matrix(X) %*% t(beta_g0)
    a <- exp(logita)/(1 + exp(logita))
    sort(apply(a,1,median)) 
    # Mean of plot means
    mean(summary(apply(a,1,mean)))
    
    # Mean predicted value on 19-April at each of the 22 sites with normal winter precipitation
    X.ppt100 <- data.frame(doy=max.doy.z,doy2=max.doy.z*max.doy.z,
                           east=g0site.locs.z[,'telem.easting'],east2=g0site.locs.z[,'telem.easting2'],
                           north=g0site.locs.z[,'telem.northing'],north2=g0site.locs.z[,'telem.northing2'],
                           ppt=(100-ppt.mn)/ppt.sd,ppt2=((100-ppt.mn)/ppt.sd)^2)
    logita.ppt100 <- matrix(1,nrow=22,ncol=1) %*% t(logita0) + as.matrix(X.ppt100) %*% t(beta_g0)
    a.ppt100 <- exp(logita.ppt100)/(1 + exp(logita.ppt100))
    sort(apply(a.ppt100,1,median))
    mean(summary(apply(a.ppt100,1,mean)))
  
  # Exploring effect of latitude
    
    northz.rng <- seq(min(cov_telem[,'telem.northing']), max(cov_telem[,'telem.northing']), length = 100)
    north.rng <- northz.rng * northing.sd + northing.mn
    northz.rng2 <- northz.rng * northz.rng
    X.n <- data.frame(north = northz.rng, north2 = northz.rng2)
    logita.n <- matrix(1,nrow=100,ncol=1) %*% t(logita0) + as.matrix(X.n) %*% t(beta_g0[,5:6])
    a.n <- exp(logita.n)/(1 + exp(logita.n))  
    a.n.med <- apply(a.n, 1, median)
    plot(a.n.med ~ north.rng)
    
#-----------------------------------------------------------------------------------------------# 
# Difference in detection distances between 2001-2003, 2004-2018
#-----------------------------------------------------------------------------------------------# 
    
  # Mean in 2001-2003 = 7.45 m
  # Mean in 2004-2018 = exp(log(mu_sigma) + beta_detect[2])
    mean(exp(log(posterior[,'mu_sigma']) + posterior[,'beta_detect[2]'])) #5.10 m
    quantile(exp(log(posterior[,'mu_sigma']) + posterior[,'beta_detect[2]']), probs = c(0.025, 0.975))

#-----------------------------------------------------------------------------------------------# 
# Marginal effects of covariates in density model
#-----------------------------------------------------------------------------------------------# 
    
  # Aspect - easting (% change for extreme values [west = -0.98 and east = 0.98])
    
    summary(surveys$aspect.e) #-0.985 to +0.980
    summary(cov_lam[,'aspect.e']) #-2.15 to +2.03
    # Density[west-facing]/Density[east-facing] = exp((-2.15-2.03)*beta)  
    aspe.diff <- min(cov_lam[,'aspect.e']) - max(cov_lam[,'aspect.e'])
    beta.east <- posterior[,'beta_lam[1]']
    east.change <- exp(aspe.diff*beta.east)
    quantile(east.change,c(0.025,0.5,0.975)) 
    
  # Aspect - northing (% change for extreme values [south = -0.98 and north = 0.98])
    
    summary(surveys$aspect.n) #-0.977 to +0.974
    summary(cov_lam[,'aspect.n']) #-2.14 to +2.11
    # Density[south-facing]/Density[north-facing] = exp((-2.14-2.11)*beta)  
    aspn.diff <- min(cov_lam[,'aspect.n']) - max(cov_lam[,'aspect.n'])
    beta.north <- posterior[,'beta_lam[2]']
    north.change <- exp(aspn.diff*beta.north)
    quantile(north.change,c(0.025,0.5,0.975))
    
  # Depth to bedrock (% increase for each x% increase in depth)
  # Note: must do % increase in bedrock because it was logged.
    
    bedL1 <- 10000
    (bedH1 <- bedL1*1.1)
    bedL1l <- log(bedL1)
    bedH1l <- log(bedH1)
    (bed10inc <- bedH1l - bedL1l)
    bedL3 <- 10000
    (bedH3 <- bedL3*1.25)
    bedL3l <- log(bedL3)
    bedH3l <- log(bedH3)
    (bed25inc <- bedH3l - bedL3l) 
    
    # Density[high] = Density[low] * exp((log(bedH) - log(bedL))*beta)
    beta.bedrock <- posterior[,'beta_lam[3]']
    # 10% increase:
    bedl.diff10 <- bed10inc/sd(log(surveys$bedrock)) #10% diff in bedrock values, on standardized scale
    bedrock.change10 <- exp(bedl.diff10*beta.bedrock)
    quantile(bedrock.change10,c(0.025,0.5,0.975))
    # 25% increase
    bedl.diff25 <- bed25inc/sd(log(surveys$bedrock)) #25% diff in bedrock values, on standardized scale
    bedrock.change25 <- exp(bedl.diff25*beta.bedrock)
    quantile(bedrock.change25,c(0.025,0.5,0.975))     

  # Wash (x% increase in density if proportion wash is halved [negative coefficient])
  # Wash was logged too (so we're evaluating a halving wash proportion)
    
    summary(surveys$wash) #0.00050 - 0.56
    logwash5 <- log(0.05)
    logwash10 <- log(0.10)
    logwash5.z <- (logwash5 - mean(log(surveys$wash)))/sd(log(surveys$wash))
    logwash10.z <- (logwash10 - mean(log(surveys$wash)))/sd(log(surveys$wash))
    beta.wash <- posterior[,'beta_lam[11]']
    wash.change <- exp((logwash5.z - logwash10.z)*beta.wash)
    quantile(wash.change,c(0.025,0.5,0.975)) 

  # Surface roughness (% increase for every 5% decrease in roughness)
    summary(surveys$rough) #1.001 - 1.190
    beta.rough <- posterior[,'beta_lam[7]']
    roughH <- 1.1
    roughL <- roughH*0.95
    rough5 <- log(roughL)-log(roughH)
    # 5% decrease
    roughl.5 <- rough5/sd(log(surveys$rough))
    rough.change5 <- exp(roughl.5*beta.rough)
    quantile(rough.change5,c(0.025,0.5,0.975)) #1.3131 (1.0210, 1.7071) so 31% increase 

  # Temp.max (identify temp.max associated with max density - quadratic term)
    summary(surveys$temp.max) #31.72 - 43.96 
    temp.max.mn <- mean(surveys$temp.max)
    temp.max.sd <- sd(surveys$temp.max)
  	#maxima = -b/2a (ax2 + bx + c)
  	b <- posterior[,'beta_lam[8]']
  	a <- posterior[,'beta_lam[9]']
  	max.temp.z <- -b/(2*a)
  	max.temp <- max.temp.z * temp.max.sd + temp.max.mn
  	quantile(max.temp,c(0.025,0.5,0.975))
