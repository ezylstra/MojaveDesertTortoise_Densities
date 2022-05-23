#===============================================================================================# 

  ## Authors: ER Zylstra, RJ Steidl, N Pope
  
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
  
    library(plyr)
    library(rstan)
    rstan_options(auto_write = TRUE)
    rstan_options(javascript = FALSE)
    library(raster)
    library(rgdal)
    library(ggplot2)

  # rm(list=ls())
  
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
    
    telem.easting2   <- telem.easting*telem.easting
    telem.northing2  <- telem.northing*telem.northing
    telem.doy2       <- telem.doy*telem.doy
    telem.precip2    <- telem.precip*telem.precip
  
#-----------------------------------------------------------------------------------------------# 
# Create indices for year, recovery units, transects
#-----------------------------------------------------------------------------------------------#  

  year <- as.integer(surveys$season-minYear + 1)   #Set minYear as year one
  surveys$recovID.ind <- as.integer(as.factor(surveys$recovID)) # CD = 1; EM = 2; NEM = 3; WM = 4  
  recovID <- surveys$recovID.ind   
  surveys$transectID.ind <- as.integer(as.factor(surveys$transectID))
  transectID <- surveys$transectID.ind
    
#-----------------------------------------------------------------------------------------------# 
# Format covariates associated with distance-sampling segments
#-----------------------------------------------------------------------------------------------#
    
  # Transform and standardize potential covariates, then save as separate objects
    # Only transforming variables when necessary to keep scaled values in a reasonable range.
    # Log transforming assumes a change in lower values has a bigger effect than a change at higher values

    precip.s <- scale(surveys$precip.s)        #Summer precip (May-Oct; dry season), 30-year normals (mm)
    precip.w <- scale(surveys$precip.w)        #Winter precip (Nov-Apr; wet season), 30-year normals (mm)
    elev <- scale(surveys$elev)                #Elevation (m)
    slope <- scale(log(surveys$slope))         #Slope (degree)
    aspect.n <- scale(surveys$aspect.n)        #Northness (1 = north; -1 = south)
    aspect.e <- scale(surveys$aspect.e)        #Eastness (1 = east; -1 = west)
    rough <- scale(log(surveys$rough))         #Average surface roughness (ratio of area to planimetric area)
    wash <- scale(log(surveys$wash))           #Proportion of grid cell classified as "wash"
    bedrock <- scale(log(surveys$bedrock))     #Depth to bedrock (cm)
    veg.p <- scale(log(surveys$veg.p))         #Perennial vegetation (NDVI value in very dry year [2019])
    veg.a <- scale(surveys$veg.a)              #Annual plant (forage) potential (comparison of MODIS values in dry, wet years)
    road.a <- scale(log(surveys$road.a))       #Distance to any road (m)
    annprecip <- scale(surveys$pptOM.perc)     #Precipitation during previous winter (Oct-Mar, % 30yr norms) -- time varying
    
    lateyrs <- 1*(surveys$season>2003)         #Indicator for all years after 2003 

  # Evaluate potential correlations:
    
    round(cor(cbind(precip.w,precip.s,
                    elev,slope,rough,aspect.n,aspect.e,
                    wash,bedrock,
                    veg.p,veg.a,
                    road.a)),3)

    round(cor(cbind(precip.w,precip.s,elev,rough,aspect.n,aspect.e,
                    wash,bedrock,veg.a,road.a)),3)  
    # With this combination, most extreme r = -0.498 (rough-bedrock)
    
  # Covariates that can potentially be used to predict availability on a given survey
    # Note: values drawn from surveys dataframe, but standardized based on mean/SD from telemetry data
    
    g0.doy <- (surveys$doy - telem.doy.mn)/telem.doy.sd
    g0.easting <- (surveys$start_e - telem.easting.mn)/telem.easting.sd
    g0.northing <- (surveys$start_n - telem.northing.mn)/telem.northing.sd
    g0.precip <- (surveys$pptOM.perc - telem.precip.mn)/telem.precip.sd
    
    g0.doy2 <- g0.doy*g0.doy
    g0.easting2 <- g0.easting*g0.easting
    g0.northing2 <- g0.northing*g0.northing
    g0.precip2 <- g0.precip*g0.precip

#-----------------------------------------------------------------------------------------------# 
# Select covariates for each "sub-model"
#-----------------------------------------------------------------------------------------------#
    
  # Choose abundance covariates at the segment-level [time-invariant]
  # cov_lam[n_segments, n_cov_lam] 

    cov_lam <- data.frame(precip.w, precip.s, elev, rough, aspect.e, aspect.n, 
                          wash, bedrock, veg.a, road.a)
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
                        precip.w,precip.s,elev,slope,rough,aspect.n,aspect.e,  
                        wash,bedrock,veg.p,veg.a,road.a)
    segCov_RUmeans <- ddply(segsc,.(recovID),summarize, 
                            precip.w=mean(precip.w),precip.s=mean(precip.s),
                            elev=mean(elev),slope=mean(slope),rough=mean(rough),
                            aspect.n=mean(aspect.n),aspect.e=mean(aspect.e),  
                            wash=mean(wash),bedrock=mean(bedrock),
                            veg.p=mean(veg.p),veg.a=mean(veg.a),
                            road.a=mean(road.a))
    
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

    ni <- 2500    #No. iterations (including warmup)
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
                  recov_l_raw=matrix(runif(n_recov*n_years,-1,1),
                                     nrow=n_recov,ncol=n_years),
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
         file="MDT_Density_2001-2020.Rdata")
  
  # Summary of model run  
    
    print(out,digits=3) 
  
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

    cexfix <- 0.9 
  
    # jpeg('Figures/Trends_2001-2020.jpg',width=165,height=228,units='mm',res=600)
    par(mfrow=c(n_recov,1),mar=c(0.5,3.0,0.5,1.0)+0.1,oma=c(2.0,0,0,0),cex=cexfix)
    # Colorado Desert
    plot(t.md~year,data=trends[trends$ru==1,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,8),bty='n',xaxs="i",yaxs='i')
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
      axis(1,at=seq(2002,2018,by=4),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0))
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
      axis(2,at=seq(0,8,by=2),labels=seq(0,8,by=2),tcl=-0.25,las=1,mgp=c(1.5,0.5,0))     
      lines(t.l~year,data=trends[trends$ru==1,],lty=2)
      lines(t.u~year,data=trends[trends$ru==1,],lty=2)
      arrows(x0=dens$year[dens$ru==1],y0=dens$D.l[dens$ru==1],x1=dens$year[dens$ru==1],y1=dens$D.u[dens$ru==1],
           col='darkgray',length=0)
      points(D.md~year,data=dens[dens$ru==1,],pch=19,cex=0.9)
      mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      text(x=2020.2,y=7.5,adj=c(1,0),'Colorado Desert',cex=0.9)
    # Eastern Mojave
    plot(t.md~year,data=trends[trends$ru==2,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,8),bty='n',xaxs="i",yaxs='i')
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
      axis(1,at=seq(2002,2018,by=4),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0))
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
      axis(2,at=seq(0,8,by=2),labels=seq(0,8,by=2),tcl=-0.25,las=1,mgp=c(1.5,0.5,0))     
      lines(t.l~year,data=trends[trends$ru==2,],lty=2)
      lines(t.u~year,data=trends[trends$ru==2,],lty=2)
      arrows(x0=dens$year[dens$ru==2],y0=dens$D.l[dens$ru==2],x1=dens$year[dens$ru==2],y1=dens$D.u[dens$ru==2],
           col='darkgray',length=0)
      points(D.md~year,data=dens[dens$ru==2,],pch=19,cex=0.9)
      mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      text(x=2020.2,y=7.5,adj=c(1,0),'Eastern Mojave',cex=0.9)
    # Northeastern Mojave
    plot(t.md~year,data=trends[trends$ru==3,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,8),bty='n',xaxs="i",yaxs='i')
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
      axis(1,at=seq(2002,2018,by=4),labels=F,tcl=-0.25,mgp=c(1.5,0.5,0))
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
      axis(2,at=seq(0,8,by=2),labels=seq(0,8,by=2),tcl=-0.25,las=1,mgp=c(1.5,0.5,0))     
      lines(t.l~year,data=trends[trends$ru==3,],lty=2)
      lines(t.u~year,data=trends[trends$ru==3,],lty=2)
      arrows(x0=dens$year[dens$ru==3],y0=dens$D.l[dens$ru==3],x1=dens$year[dens$ru==3],y1=dens$D.u[dens$ru==3],
           col='darkgray',length=0)
      points(D.md~year,data=dens[dens$ru==3,],pch=19,cex=0.9)
      mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      text(x=2020.2,y=7.5,adj=c(1,0),'Northeastern Mojave',cex=0.9)
    # Western Mojave
    plot(t.md~year,data=trends[trends$ru==4,],type='l',xaxt='n',yaxt='n',xlab='',ylab='',
         xlim=c(2000.5,2020.5),ylim=c(-0.3,8),bty='n',xaxs="i",yaxs='i')
      usr <- par('usr')  #these are plotting limits (incl extra bit)
      axis(1,at=c(usr[1],usr[2]),tck=F,labels=F)
      axis(1,at=seq(2002,2018,by=4),labels=seq(2002,2018,by=4),tcl=-0.25,mgp=c(1.5,0.5,0))
      axis(2,at=c(usr[3],usr[4]),tck=F,labels=F)
      axis(2,at=seq(0,8,by=2),labels=seq(0,8,by=2),tcl=-0.25,las=1,mgp=c(1.5,0.5,0))     
      lines(t.l~year,data=trends[trends$ru==4,],lty=2)
      lines(t.u~year,data=trends[trends$ru==4,],lty=2)
      arrows(x0=dens$year[dens$ru==4],y0=dens$D.l[dens$ru==4],x1=dens$year[dens$ru==4],y1=dens$D.u[dens$ru==4],
           col='darkgray',length=0)
      points(D.md~year,data=dens[dens$ru==4,],pch=19,cex=0.9)
      mtext('Density (tortoises/sq.km)',side=2,las=0,line=2.0,cex=cexfix)
      mtext('Year',side=1,line=1.5,cex=cexfix)
      text(x=2020.2,y=7.5,adj=c(1,0),'Western Mojave',cex=0.9)
    # dev.off()
  
#-----------------------------------------------------------------------------------------------# 
# Calculating marginal effects of covariates on predicted availability 
#-----------------------------------------------------------------------------------------------# 
  
  easting.mn <- 653940.5
  easting.sd <- 71116.53
  doy.mn <- 114.6716
  doy.sd <- 18.67406
  ppt.mn <- 96.40828
  ppt.sd <- 62.21435

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
    
    # Mean predicted value on 20-April at each of the 22 sites with maximum winter precipitation
    X <- data.frame(doy=max.doy.z,doy2=max.doy.z*max.doy.z,
                    east=g0site.locs.z[,'telem.easting'],east2=g0site.locs.z[,'telem.easting2'],
                    north=g0site.locs.z[,'telem.northing'],north2=g0site.locs.z[,'telem.northing2'],
                    ppt=max.ppt.z,ppt2=max.ppt.z*max.ppt.z)
    logita <- matrix(1,nrow=22,ncol=1) %*% t(logita0) + as.matrix(X) %*% t(beta_g0)
    a <- exp(logita)/(1 + exp(logita))
    sort(apply(a,1,median)) 
    # Mean of plot means
    mean(summary(apply(a,1,mean)))
    
    # Mean predicted value on 20-April at each of the 22 sites with normal winter precipitation
    X.ppt100 <- data.frame(doy=max.doy.z,doy2=max.doy.z*max.doy.z,
                           east=g0site.locs.z[,'telem.easting'],east2=g0site.locs.z[,'telem.easting2'],
                           north=g0site.locs.z[,'telem.northing'],north2=g0site.locs.z[,'telem.northing2'],
                           ppt=(100-ppt.mn)/ppt.sd,ppt2=((100-ppt.mn)/ppt.sd)^2)
    logita.ppt100 <- matrix(1,nrow=22,ncol=1) %*% t(logita0) + as.matrix(X.ppt100) %*% t(beta_g0)
    a.ppt100 <- exp(logita.ppt100)/(1 + exp(logita.ppt100))
    sort(apply(a.ppt100,1,median))
    # Mean of plot means
    mean(summary(apply(a.ppt100,1,mean)))
  
#-----------------------------------------------------------------------------------------------# 
# Difference in detection distances between 2001-2003, 2004-2018
#-----------------------------------------------------------------------------------------------# 
    
  # Mean in 2001-2003 = 7.42 m
  # Mean in 2004-2018 = exp(log(mu_sigma) + beta_detect[2])
    mean(exp(log(posterior[,'mu_sigma']) + posterior[,'beta_detect[2]'])) #5.08 m

#-----------------------------------------------------------------------------------------------# 
# Marginal effects of covariates in density model
#-----------------------------------------------------------------------------------------------# 
    
  # Aspect - easting (% change for extreme values [west = -0.98 and east = 0.98])
    
    summary(surveys$aspect.e) #-0.985 to +0.980
    summary(cov_lam[,'aspect.e']) #-2.15 to +2.03
    #Density[west-facing] = Density[east-facing] * exp((-2.15-2.03)*beta  
    aspe.diff <- min(cov_lam[,'aspect.e']) - max(cov_lam[,'aspect.e'])
    beta.east <- posterior[,'beta_lam[1]']
    east.change <- exp(aspe.diff*beta.east)
    quantile(east.change,c(0.025,0.5,0.975)) #1.36 (1.18, 1.57), so 36% increase
    
  # Aspect - northing (% change for extreme values [south = -0.98 and north = 0.98])
    
    summary(surveys$aspect.n) #-0.977 to +0.974
    summary(cov_lam[,'aspect.n']) #-2.14 to +2.11
    #Density[south-facing] = Density[north-facing] * exp((-2.14-2.11)*beta  
    aspn.diff <- min(cov_lam[,'aspect.n']) - max(cov_lam[,'aspect.n'])
    beta.north <- posterior[,'beta_lam[2]']
    north.change <- exp(aspn.diff*beta.north)
    quantile(north.change,c(0.025,0.5,0.975)) #1.16 (1.01, 1.34), so 16% increase    
    
  # Depth to bedrock (% increase for each x% increase in depth)
  # Note: must calculate % increase in bedrock because it was logged.
    
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
    bedl.diff10 <- bed10inc/sd(log(surveys$bedrock)) # 10% diff in bedrock values, on standardized scale
    bedrock.change10 <- exp(bedl.diff10*beta.bedrock)
    quantile(bedrock.change10,c(0.025,0.5,0.975)) # 1.005 (0.999, 1.012), so 0.5% increase   
    # 25% increase
    bedl.diff25 <- bed25inc/sd(log(surveys$bedrock)) # 25% diff in bedrock values, on standardized scale
    bedrock.change25 <- exp(bedl.diff25*beta.bedrock)
    quantile(bedrock.change25,c(0.025,0.5,0.975)) # 1.0121 (0.9976, 1.0271), so 1.2% increase       

  # Wash (x% increase in density between wash = 0.05 and wash = 0.10)
  # Wash was logged too (so we're evaluating a doubling in wash proportion)
    
    summary(surveys$wash)
    logwash5 <- log(0.05)
    logwash10 <- log(0.10)
    logwash5.z <- (logwash5 - mean(log(surveys$wash)))/sd(log(surveys$wash))
    logwash10.z <- (logwash10 - mean(log(surveys$wash)))/sd(log(surveys$wash))
    beta.wash <- posterior[,'beta_lam[10]']
    wash.change <- exp((logwash5.z - logwash10.z)*beta.wash)
    quantile(wash.change,c(0.025,0.5,0.975)) #1.0226 (0.9944, 1.0507), so 2.3% increase
    
  # Roads (% increase for every 50% increase in distance)
    
    summary(surveys$road.a)
    beta.road <- posterior[,'beta_lam[7]']
    roadL <- 1000
    roadH <- roadL*1.50
    road50 <- log(roadH)-log(roadL)
    # 50% increase
    roadl.50 <- road50/sd(log(surveys$road.a))
    road.change50 <- exp(roadl.50*beta.road)
    quantile(road.change50,c(0.025,0.5,0.975)) #1.0095 (0.9995, 1.0202), so 0.95% increase         

  # Surface roughness (% increase for every 5% decrease in roughness)
    
    summary(surveys$rough) #1.001 - 1.190
    beta.rough <- posterior[,'beta_lam[8]']
    roughH <- 1.1
    roughL <- roughH*0.95
    rough5 <- log(roughL)-log(roughH)
    # 5% decrease
    roughl.5 <- rough5/sd(log(surveys$rough))
    rough.change5 <- exp(roughl.5*beta.rough)
    quantile(rough.change5,c(0.025,0.5,0.975)) #1.3477 (1.0406, 1.7519) so 35% increase 
  
#-----------------------------------------------------------------------------------------------# 
# Map with survey locations, telemetry plots, prediction area
#-----------------------------------------------------------------------------------------------# 

  # Load shapefiles and convert to dataframes for ggplot
    
    # One covariate raster (to get projection)
    
      east <- raster('Covariates/avg_eastness.tif')
      
    # Recovery units
      
      rus <- readOGR(dsn='Covariates/Revised Recovery Units',layer='2011RecoveryUnits')
      rus <- spTransform(rus,crs(east))
      rus@data$id <- rownames(rus@data)
      rus_data <- fortify(rus, unit = "id")
      rus_df <- merge(rus_data, rus@data, by = "id")
      rus_albers <- spTransform(rus, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 
                                         +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
    
    # TCAs
      
      tcas <- readOGR(dsn='Covariates/TCAs',layer='All_Strata')
      tcas <- spTransform(tcas,crs(east))
      tcas_noUVR <- subset(tcas, recovery_u != "Upper Virgin River")
      tcas_noUVR@data$id <- rownames(tcas_noUVR@data)
      tcas_data <- fortify(tcas_noUVR, unit = "id")
      tcas_df <- merge(tcas_data, tcas_noUVR@data, by = "id")
    
    # State boundaries
      
      states <- readOGR(dsn='Covariates/US_states_GIS',layer='cb_2017_us_state_500k')
      states <- spTransform(states,crs(east))
      cstates <- subset(states, !STUSPS %in% c("MP", "AS", "PR", "VI", "HI" , "AK", "GU"))
      cstates@data$id <- rownames(cstates@data)
      cstates_data <- fortify(cstates, unit = "id")
      cstates_df <- merge(cstates_data, cstates@data, by = "id")
      cstates_albers <- spTransform(cstates, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 
                                                 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

    # Prediction Area (includes all 1 sq km grid cells in 4 RUs EXCEPT those cells where one or
      # more covariate values is > 10% outside the range of values at survey locations)

      predarea <- readOGR(dsn='Covariates/PredictionArea_Full',layer='PredictionArea_Full')
      predarea <- spTransform(predarea,crs(east))  
      predarea@data$id <- rownames(predarea@data)
      predarea_data <- fortify(predarea, unit = "id")
      predarea_df <- merge(predarea_data, predarea@data, by = "id")
        
    # Highways
      
      roads.nv <- readOGR(dsn='Covariates/Roads_NV',layer='tl_2021_32_prisecroads')
      roads.nv <- spTransform(roads.nv,crs(east))
      roads.ca <- readOGR(dsn='Covariates/Roads_CA',layer='tl_2021_06_prisecroads')
      roads.ca <- spTransform(roads.ca,crs(east))
      roads.az <- readOGR(dsn='Covariates/Roads_AZ',layer='tl_2021_04_prisecroads')
      roads.az <- spTransform(roads.az,crs(east))
      roads.ut <- readOGR(dsn='Covariates/Roads_UT',layer='tl_2021_49_prisecroads')
      roads.ut <- spTransform(roads.ut,crs(east))

      highways.nv <- subset(roads.nv,MTFCC=='S1100')
      highways.ca <- subset(roads.ca,MTFCC=='S1100')
      highways.az <- subset(roads.az,MTFCC=='S1100')
      highways.ut <- subset(roads.ut,MTFCC=='S1100')
      
      highways.az@data$id <- rownames(highways.az@data)
      highways.az_data <- fortify(highways.az, unit = "id")
      highways.az_df <- merge(highways.az_data, highways.az@data, by = "id")
      highways.nv@data$id <- rownames(highways.nv@data)
      highways.nv_data <- fortify(highways.nv, unit = "id")
      highways.nv_df <- merge(highways.nv_data, highways.nv@data, by = "id")
      highways.ut@data$id <- rownames(highways.ut@data)
      highways.ut_data <- fortify(highways.ut, unit = "id")
      highways.ut_df <- merge(highways.ut_data, highways.ut@data, by = "id")
      highways.ca@data$id <- rownames(highways.ca@data)
      highways.ca_data <- fortify(highways.ca, unit = "id")
      highways.ca_df <- merge(highways.ca_data, highways.ca@data, by = "id")
      
    # Extract survey locations (segment midpoints)
    
      midlocs <- surveys[,c('mid_easting','mid_northing')]

    # Extract telemetry plot locations
      
      plotlocs <- unique(g0obs[,c('site','east','north')])

  options(scipen =1000000)
      
  gSA <- ggplot() +
          geom_polygon(data = predarea_df, aes(x = long, y = lat, group = group), 
                       fill = "gray85", color = NA, size = NA) + 
          geom_point(data = midlocs, aes(x = mid_easting, y = mid_northing), 
                     color = "steelblue3", size = 0.3) +          
          geom_path(data = highways.ca_df, aes(x = long, y = lat, group = group),
                    color = "gray50", size = 0.15) +
          geom_path(data = highways.ut_df, aes(x = long, y = lat, group = group),
                    color = "gray50", size = 0.15) +
          geom_path(data = highways.nv_df, aes(x = long, y = lat, group = group),
                    color = "gray50", size = 0.15) +
          geom_path(data = highways.az_df, aes(x = long, y = lat, group = group),
                    color = "gray50", size = 0.15) +
          geom_path(data = cstates_df, aes(x = long, y = lat, group = group),
                    color = "gray50", size = 0.7) +
          geom_path(data = rus_df, aes(x = long, y = lat, group = group), 
                    color = "black", size = 1) +      
          geom_path(data = tcas_df, aes(x = long, y = lat, group = group), 
                    color = "black", size = 0.5) + 
          geom_point(data = plotlocs, aes(x = east, y = north), 
                     fill = "yellow", color = "black", pch = 23, size = 2) + 
          annotate("text", x = 470000, y = 3970000, label = "WM", hjust = 0.5, 
                   vjust = 0.5, size = 12/.pt, fontface = 2) +
          annotate("text", x = 580000, y = 4080000, label = "EM", hjust = 0.5, 
                   vjust = 0.5, size = 12/.pt, fontface = 2) +
          annotate("text", x = 700000, y = 4138000, label = "NM", hjust = 0.5, 
                   vjust = 0.5, size = 12/.pt, fontface = 2) +
          annotate("text", x = 700000, y = 3752000, label = "CD", hjust = 0.5, 
                   vjust = 0.5, size = 12/.pt, fontface = 2) +
          annotate("text", x = 805000, y = 4120000, label = "UVR", hjust = 0.5, 
                   vjust = 0.5, size = 12/.pt, fontface = 2) +
          theme(axis.title = element_blank(),
                axis.text.y = element_text(angle = 90, hjust = 0.5),
                panel.background = element_blank(),
                panel.border = element_rect(color = 'black', fill = NA),
                plot.margin = grid::unit(c(0.1,0,0,0), "in")) + 
          coord_fixed(ratio = 1, xlim=c(336500, 825000), ylim=c(3631000, 4155000))
  
  ginset <- ggplot() + 
              geom_polygon(data = cstates_albers, aes(x = long, y = lat, group = group), 
                                     fill = NA, color = "gray30", size = 0.2) + 
              geom_polygon(data = rus_albers, aes(x = long, y = lat, group = group), fill = "black") + 
              theme(axis.title = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(color = 'black', fill = NA)) +
              coord_fixed(ratio = 1)
  
  gSA + annotation_custom(grob=ggplotGrob(ginset), 
                          xmin = 310000, xmax = 510000, ymin = 3600000, ymax = 3760000) 
  
  ggsave("Figures/StudyArea_wInset.jpg",
         device="jpeg",
         height = 6.8,
         width = 6.5,
         units = "in",
         dpi = 600)
      
