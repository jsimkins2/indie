##' Create statistical models to predict subdaily meteorology from daily means
##' 
##' @title gen_subdaily_models
##' @param outfolder - directory where models will be stored *** storage required varies by size of training dataset, but prepare for >100 GB
##' @param in.path - filepath to training dataset created by CF2traindata.R
##' @param in.prefix - prefix of train_data, i.e. if file is US-WCr_train_data, prefix is "US-WCr"
##' @param tdf_filepath - temporal_downscale_functions.R filepath that can be sourced i.e. "~/scripts/temporal_downscale_functions.R"
##' @param n.beta - number of betas to save from linear regression model
##' @param resids - logical stating whether to pass on residual data or not
##' @param parallel - logical stating whether to run temporal_downscale_functions.R in parallel 
##' @param n.cores - deals with parallelization
##' @param day.window - integer specifying number of days around a particular day you want to use data from for that 
##'                     specific hours coefficients
##' @param overwrite
##' @param verbose

##' @author Christy Rollinson, James Simkins

gen_subdaily_models <- function(outfolder, in.path, in.prefix, tdf_filepath, n.beta,
                                resids=F, parallel=F, n.cores=NULL, day.window,
                                overwrite = TRUE, verbose = FALSE){
  # ------------------------------------------
  # 0 Load Libraries, set up file directories
  # ------------------------------------------
  # Script to prototype temporal downscaling

  # ------------------------------------------
  # 1. Load and format training data
  # ------------------------------------------
  # ----------
  # 1.0 Read data & Make time stamps
  # ----------
    # Load the data
  train_path = file.path(in.path, paste0(in.prefix, "_train_data"))
  dat.train <- read.csv(train_path)
  
  dat.train$date <- strptime(paste(dat.train$year, dat.train$doy+1, dat.train$hour, sep="-"), "%Y-%j-%H", tz="GMT")
  dat.train$time.hr <- as.numeric(difftime(dat.train$date, paste0((min(dat.train$year)-1),"-12-31 23:00:00"), tz="GMT", units="hour"))
  dat.train$time.day <- as.numeric(difftime(dat.train$date, paste0((min(dat.train$year)-1),"-12-31 23:00:00"), tz="GMT", units="day"))-1/24
  dat.train$time.day2 <- as.integer(dat.train$time.day+1/(48*2))+1 # Offset by half a time step to get time stamps to line up
  dat.train <- dat.train[order(dat.train$time.hr),]
  

  # ----------
  # 1.1 Coming up with the daily means that are what we can use as predictors
  # ----------
  
  train.day <- aggregate(dat.train[,c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")], 
                         by=dat.train[,c("year", "doy")],
                         FUN=mean)
  names(train.day)[3:9] <- c("tmean.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")
  train.day$tmax.day <- aggregate(dat.train[,c("tair")], by=dat.train[,c("year", "doy")], FUN=max)$x
  train.day$tmin.day <- aggregate(dat.train[,c("tair")], by=dat.train[,c("year", "doy")], FUN=min)$x
  summary(train.day)
  
  
  dat.train <- merge(dat.train[,], train.day, all.x=T, all.y=T)
  summary(dat.train)

# ----------
# 1.2 Setting up a 1-hour lag -- smooth transitions at midnight
# NOTE: because we're filtering from the present back through the past, -1 will associate the closest 
#       hour that we've already done (midnight) with the day we're currently working on
# ----------
  vars.hour <- c("tair","precipf", "swdown", "lwdown", "press", "qair", "wind")
  vars.lag <- c("lag.tair", "lag.precipf", "lag.swdown", "lag.lwdown", "lag.press", "lag.qair", "lag.wind") 
  lag.day <- dat.train[dat.train$hour==23,c("year", "doy", "time.day2", vars.hour)]
  names(lag.day)[4:10] <- vars.lag

  lag.day <- aggregate(lag.day[,vars.lag],
                       by=lag.day[,c("year", "doy", "time.day2")],
                       FUN=mean)
  lag.day$lag.tmin <- aggregate(dat.train[,c("tair")],  
                                by=dat.train[,c("year", "doy", "time.day2")],
                                FUN=min)[,"x"] # Add in a lag for the next day's min temp
  lag.day$lag.tmax <- aggregate(dat.train[,c("tair")],  
                                by=dat.train[,c("year", "doy", "time.day2")],
                                FUN=max)[,"x"] # Add in a lag for the next day's min temp
  lag.day$time.day2 <- lag.day$time.day2+1 # +1 for forward filtering downscale

  
  dat.train <- merge(dat.train, lag.day[,c("time.day2", vars.lag, "lag.tmin", "lag.tmax")], all.x=T)

# ----------
# 1.3 Setting up a variable to 'preview' the next day's mean to help get smoother transitions
# NOTE: because we're filtering from the present back through the past, +1 will associate 
#       the mean for the next day we're going to model with the one we're currently working on
# ----------
  vars.day <- c("tmean.day", "tmax.day", "tmean.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")
  vars.next <- c("next.tmean", "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind") 
  
  next.day <- dat.train[c("year", "doy", "time.day2", vars.day)]
  names(next.day)[4:12] <- vars.next
  next.day <- aggregate(next.day[,vars.next],
                        by=next.day[,c("year", "doy", "time.day2")],
                        FUN=mean)
  next.day$time.day2 <- next.day$time.day2-1  
  
  dat.train <- merge(dat.train, next.day[,c("time.day2", vars.next)], all.x=T)

# ----------
# 1.4 calculate tmin & tmax as departure from mean; order data
# ----------
  # Lookign at max & min as departure from mean
  dat.train$max.dep <- dat.train$tmax.day - dat.train$tmean.day
  dat.train$min.dep <- dat.train$tmin.day - dat.train$tmean.day
  summary(dat.train)
  
  # Order the data just to help with my sanity when visualizing
  dat.train <- dat.train[order(dat.train$time.hr),]
  summary(dat.train)
  
  # ------------------------------------------
  # 2 Train the models for each variable and save them to be read in as needed
  # ------------------------------------------
  source(tdf_filepath)
  
  # ---------
  # 2.1 Generating all the daily models, save the output as .Rdata files, then clear memory
  # Note: Could save Betas as .nc files that we pull from as needed to save memory; but for now just leaving it in the .Rdata file for eas
  # Note: To avoid propogating too much wonkiness in hourly data, any co-variates are at the daily level
  # Note: If mod.precipf.doy doesn't run, try increasing the day.window for this variable. The lack of non-zero
  #       values makes it difficult for the linear regression model to calculate coefficients sometimes
  # ---------
  
  mod.tair.doy    <- model.tair   (dat.train=dat.train[,], resids=resids, parallel=parallel,path.out=paste0(outfolder,in.prefix, "/model.tair"), n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.tair.doy)
  
  mod.precipf.doy <- model.precipf(dat.train=dat.train[,], resids=resids, parallel=parallel,path.out=paste0(outfolder,in.prefix, "/model.precipf"), n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.precipf.doy)
  
  mod.swdown.doy  <- model.swdown (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/model.swdown"),n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.swdown.doy)
  
  mod.lwdown.doy  <- model.lwdown (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/model.lwdown"),n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.lwdown.doy)
  
  mod.press.doy   <- model.press  (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/model.press"),n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.press.doy)
  
  mod.qair.doy    <- model.qair   (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/model.qair"),n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.qair.doy)
  
  mod.wind.doy    <- model.wind   (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/model.wind"),n.cores=n.cores, n.beta=n.beta, day.window=5)
  rm(mod.wind.doy)
}
