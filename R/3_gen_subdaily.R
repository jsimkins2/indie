##' Create statistical models to predict subdaily meteorology from daily means
##' 
##' @title gen_subdaily_models
##' @param outfolder - directory where models will be stored *** storage required varies by size of training dataset, but prepare for >100 GB
##' @param dat.train_file - train_data file
##' @param tdf_file - temporal_downscale_functions.R filepath that can be sourced i.e. "~/scripts/temporal_downscale_functions.R"
##' @param n.beta - number of betas to save from linear regression model
##' @param resids - logical stating whether to pass on residual data or not
##' @param parallel - logical stating whether to run temporal_downscale_functions.R in parallel 
##' @param n.cores - deals with parallelization
##' @param day.window - integer specifying number of days around a particular day you want to use data from for that 
##'                     specific hours coefficients
##' @param overwrite
##' @param verbose

##' @author Christy Rollinson, James Simkins

gen_subdaily_models <- function(outfolder, dat.train_file, tdf_file, n.beta,
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
  
  vars.info <- data.frame(CF.name = c("date", "year", "doy","hour", "air_temperature", "precipitation_flux", "air_temperature_max", "air_temperature_min",
                                      "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air","air_pressure",
                                      "specific_humidity", "eastward_wind", "northward_wind", "wind_speed"))
  dat.train <- list()
  tem <- ncdf4::nc_open(dat.train_file)
  dim <- tem$dim
  for (j in seq_along(vars.info$CF.name)) {
    if (exists(as.character(vars.info$CF.name[j]), tem$var) == FALSE) {
      dat.train[[j]] = NA
    } else {
      dat.train[[j]] = ncdf4::ncvar_get(tem, as.character(vars.info$CF.name[j]))
    }
  }
  names(dat.train)<- vars.info$CF.name 
  dat.train <- data.frame(dat.train)
  
  dat.train$date <- strptime(paste(dat.train$year, dat.train$doy+1, dat.train$hour, sep="-"), "%Y-%j-%H", tz="GMT")
  dat.train$time.hr <- as.numeric(difftime(dat.train$date, paste0((min(dat.train$year)-1),"-12-31 23:00:00"), tz="GMT", units="hour"))
  dat.train$time.day <- as.numeric(difftime(dat.train$date, paste0((min(dat.train$year)-1),"-12-31 23:00:00"), tz="GMT", units="day"))-1/24
  dat.train$time.day2 <- as.integer(dat.train$time.day+1/(48*2))+1 # Offset by half a time step to get time stamps to line up
  dat.train <- dat.train[order(dat.train$time.hr),]
  

  # ----------
  # 1.1 Coming up with the daily means that are what we can use as predictors
  # ----------
  
  train.day <- aggregate(dat.train[,c("air_temperature", "precipitation_flux", "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air", "air_pressure", "specific_humidity", "wind_speed")], 
                         by=dat.train[,c("year", "doy")],
                         FUN=mean)
  names(train.day)[3:9] <- c("air_temperature_mean.day", "precipitation_flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day")
  train.day$air_temperature_max.day <- aggregate(dat.train[,c("air_temperature")], by=dat.train[,c("year", "doy")], FUN=max)$x
  train.day$air_temperature_min.day <- aggregate(dat.train[,c("air_temperature")], by=dat.train[,c("year", "doy")], FUN=min)$x
  summary(train.day)
  
  
  dat.train <- merge(dat.train[,], train.day, all.x=T, all.y=T)
  summary(dat.train)

# ----------
# 1.2 Setting up a 1-hour lag -- smooth transitions at midnight
# NOTE: because we're filtering from the present back through the past, -1 will associate the closest 
#       hour that we've already done (midnight) with the day we're currently working on
# ----------
  vars.hour <- c("air_temperature","precipitation_flux", "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air", "air_pressure", "specific_humidity", "wind_speed")
  vars.lag <- c("lag.air_temperature", "lag.precipitation_flux", "lag.surface_downwelling_shortwave_flux_in_air", "lag.surface_downwelling_longwave_flux_in_air", "lag.air_pressure", "lag.specific_humidity", "lag.wind_speed") 
  lag.day <- dat.train[dat.train$hour==23,c("year", "doy", "time.day2", vars.hour)]
  names(lag.day)[4:10] <- vars.lag

  lag.day <- aggregate(lag.day[,vars.lag],
                       by=lag.day[,c("year", "doy", "time.day2")],
                       FUN=mean)
  lag.day$lag.air_temperature_min <- aggregate(dat.train[,c("air_temperature")],  
                                by=dat.train[,c("year", "doy", "time.day2")],
                                FUN=min)[,"x"] # Add in a lag for the next day's min temp
  lag.day$lag.air_temperature_max <- aggregate(dat.train[,c("air_temperature")],  
                                by=dat.train[,c("year", "doy", "time.day2")],
                                FUN=max)[,"x"] # Add in a lag for the next day's min temp
  lag.day$time.day2 <- lag.day$time.day2+1 # +1 for forward filtering downscale

  
  dat.train <- merge(dat.train, lag.day[,c("time.day2", vars.lag, "lag.air_temperature_min", "lag.air_temperature_max")], all.x=T)

# ----------
# 1.3 Setting up a variable to 'preview' the next day's mean to help get smoother transitions
# NOTE: because we're filtering from the present back through the past, +1 will associate 
#       the mean for the next day we're going to model with the one we're currently working on
# ----------
  vars.day <- c("air_temperature_mean.day", "air_temperature_max.day", "air_temperature_mean.day", "precipitation_flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day")
  vars.next <- c("next.air_temperature_mean", "next.air_temperature_max", "next.air_temperature_min", "next.precipitation_flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed") 
  
  next.day <- dat.train[c("year", "doy", "time.day2", vars.day)]
  names(next.day)[4:12] <- vars.next
  next.day <- aggregate(next.day[,vars.next],
                        by=next.day[,c("year", "doy", "time.day2")],
                        FUN=mean)
  next.day$time.day2 <- next.day$time.day2-1  
  
  dat.train <- merge(dat.train, next.day[,c("time.day2", vars.next)], all.x=T)

# ----------
# 1.4 calculate air_temperature_min & air_temperature_max as departure from mean; order data
# ----------
  # Lookign at max & min as departure from mean
  dat.train$max.dep <- dat.train$air_temperature_max.day - dat.train$air_temperature_mean.day
  dat.train$min.dep <- dat.train$air_temperature_min.day - dat.train$air_temperature_mean.day
  summary(dat.train)
  
  # Order the data just to help with my sanity when visualizing
  dat.train <- dat.train[order(dat.train$time.hr),]
  summary(dat.train)
  
  # ------------------------------------------
  # 2 Train the models for each variable and save them to be read in as needed
  # ------------------------------------------
  source(tdf_file)

  
  # ---------
  # 2.1 Generating all the daily models, save the output as .Rdata files, then clear memory
  # Note: Could save Betas as .nc files that we pull from as needed to save memory; but for now just leaving it in the .Rdata file for eas
  # Note: To avoid propogating too much wonkiness in hourly data, any co-variates are at the daily level
  # Note: If mod.precipitation_flux.doy doesn't run, try increasing the day.window for this variable. The lack of non-zero
  #       values makes it difficult for the linear regression model to calculate coefficients sometimes
  # ---------
  
  mod.air_temperature.doy    <- model.air_temperature   (dat.train=dat.train[,], resids=resids, parallel=parallel,path.out=paste0(outfolder,in.prefix, "/air_temperature"), n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.air_temperature.doy)

  mod.precipitation_flux.doy <- model.precipitation_flux(dat.train=dat.train[,], resids=resids, parallel=parallel,path.out=paste0(outfolder,in.prefix, "/precipitation_flux"), n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.precipitation_flux.doy)

  mod.surface_downwelling_shortwave_flux_in_air.doy  <- model.surface_downwelling_shortwave_flux_in_air (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/surface_downwelling_shortwave_flux_in_air"),n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.surface_downwelling_shortwave_flux_in_air.doy)

  mod.surface_downwelling_longwave_flux_in_air.doy  <- model.surface_downwelling_longwave_flux_in_air (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/surface_downwelling_longwave_flux_in_air"),n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.surface_downwelling_longwave_flux_in_air.doy)

  mod.air_pressure.doy   <- model.air_pressure  (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/air_pressure"),n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.air_pressure.doy)
  
  mod.specific_humidity.doy    <- model.specific_humidity   (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/specific_humidity"),n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.specific_humidity.doy)
  
  mod.wind_speed.doy    <- model.wind_speed   (dat.train=dat.train[,], resids=resids, parallel=parallel, path.out=paste0(outfolder,in.prefix, "/wind_speed"),n.cores=n.cores, n.beta=n.beta, day.window=day.window)
  rm(mod.wind_speed.doy)
}
