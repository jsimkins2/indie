##' Create statistical models to predict subdaily meteorology from daily means
##' 
##' @title gen_subdaily_models
##' @param outfolder - directory where output file will be stored
##' @param in.path - path to model dataset you wish to temporally downscale
##' @param in.prefix - prefix of model dataset, i.e. if file is GFDL.CM3.rcp45.r1i1p1.2006 the prefix is "GFDL.CM3.rcp45.r1i1p1"
##' @param dat.train_file - location of train_data file
##' @param lm.models.base - path to linear regression model folder from 3_gen_subdaily
##' @param start_date
##' @param end_date
##' @param tdf_file - temporal_downscale_functions.R filepath that can be sourced i.e. "~/scripts/temporal_downscale_functions.R"
##' @param td_file - temporal_downscale.R filepath that can be sourced i.e. "~/scripts/temporal_downscale.R"
##' @param site.name - name of location i.e. Ameriflux Site as character string
##' @param lat.in - latitude as numeric
##' @param lon.in - longitude as numeric
##' @param ens.hr - integer selecting number of hourly ensemble members 
##' @param n.day - 1 # Number of daily ensemble members to process, must be 1 because not passing in daily ensemble
##' @param cores.max - 12
##' @param resids - logical stating whether to pass on residual data or not
##' @param parallel - logical stating whether to run temporal_downscale_functions.R in parallel 
##' @param n.cores - deals with parallelization
##' @param day.window - integer specifying number of days around a particular day you want to use data from for that 
##'                     specific hours coefficients
##' @param overwrite
##' @param verbose

##' @author Christy Rollinson, James Simkins
##' 

pred_subdaily_met <- function(outfolder, in.path, in.prefix, lm.models.base, dat.train_file, 
                              start_date, end_date, tdf_file, td_file, site.name, lat.in, lon.in,  
                              cores.max = 12, ens.hr = 3, n.day = 1, resids = FALSE, 
                              parallel = FALSE, n.cores=NULL) {
  
  years = seq(lubridate::year(start_date), lubridate::year(end_date))

  # Load the scripts that do all the heavy lifting
  source(tdf_file)
  source(td_file)

  # Load the training dataset
  dat.train <- read.csv(dat.train_file)
  
  df.hour <- data.frame(hour=unique(dat.train$hour)) # match this to whatever your "hourly" timestep is
  
  # Set up the appropriate seed
  set.seed(0017)
  seed.vec <- sample.int(1e6, size=500, replace=F)
  
  # Defining variable names, longname & units
  vars.info <- data.frame(name    =c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind"),
                          name.cf = c("air_temperature", 
                                      "precipitation_flux",
                                      "surface_downwelling_shortwave_flux_in_air",
                                      "surface_downwelling_longwave_flux_in_air",
                                      "air_pressure",
                                      "specific_humidity",
                                      "wind"
                          ),
                          longname=c("2 meter mean air temperature", 
                                     "cumulative precipitation (water equivalent)",
                                     "incident (downwelling) showtwave radiation",
                                     "incident (downwelling) longwave radiation",
                                     'Pressure at the surface',
                                     'Specific humidity measured at the lowest level of the atmosphere',
                                     'Wind speed' 
                          ),
                          units= c("K", "kg m-2 s-1", "W m-2", "W m-2", "Pa", "kg kg-1", "m s-1")
  )
  # ----------------------------------
  for (y in years){

    path.gcm <- file.path(in.path, paste0(in.prefix,".", y, ".nc"))
    path.out <- file.path(in.path, paste0("hr"))
    if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
    
    # -----------------------------------
    # 1. Format output so all ensemble members can be run at once
    # NOTE: Need to start with the last and work to the first
    # -----------------------------------
    # Figure out which daily ensembles we should pull from
    day.dirs <- dir(path.gcm,  paste0(site.name, "_", in.prefix, "_day_")) #directory of input data aka maca
    ens.day.all <- substr(day.dirs, nchar(day.dirs)-2, nchar(day.dirs))
    
    hrs.dir <- dir(file.path(path.out),  paste0(site.name, "_", in.prefix, "_1hr_")) #this is the output path
    ens.day.done <- substr(hrs.dir, nchar(hrs.dir)-6, nchar(hrs.dir)-4)
    
    ens.list <- ens.day.all[!(ens.day.all %in% ens.day.done)] #list of files so aka multiple years
    
    seed <- seed.vec[length(hrs.dir)+1] # This makes sure that if we add ensemble members, it gets a new, but reproducible seed
    set.seed(seed)
    
    
    
    # Initialize the lags
    lags.init <- list() # Need to initialize lags first outside of the loop and then they'll get updated internally
    # by doing this, will they get updated in temp_dwnsc.??? 
  
      
    nc.now <- ncdf4::nc_open(path.gcm)
    # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
    nc.time <- ncdf4::ncvar_get(nc.now, "time")
    tmax.init    <- ncdf4::ncvar_get(nc.now, "air_temperature_max")[length(nc.time)]
    tmin.init    <- ncdf4::ncvar_get(nc.now, "air_temperature_min")[length(nc.time)]
    precipf.init <- ncdf4::ncvar_get(nc.now, "precipitation_flux")[length(nc.time)]
    swdown.init  <- ncdf4::ncvar_get(nc.now, "surface_downwelling_shortwave_flux_in_air")[length(nc.time)]
    lwdown.init  <- ncdf4::ncvar_get(nc.now, "surface_downwelling_longwave_flux_in_air")[length(nc.time)]
    press.init   <- ncdf4::ncvar_get(nc.now, "air_pressure")[length(nc.time)]
    qair.init    <- ncdf4::ncvar_get(nc.now, "specific_humidity")[length(nc.time)]
    nwind.init    <- ncdf4::ncvar_get(nc.now, "northward_wind")[length(nc.time)]
    ewind.init <- ncdf4::ncvar_get(nc.now, "eastward_wind")[length(nc.time)]
    wind.init <- sqrt(nwind.init^2 + ewind.init^2)
    ncdf4::nc_close(nc.now)

    lags.init[["tair"   ]] <- data.frame(array(mean(c(tmax.init, tmin.init)), dim=c(1, ens.hr)))
    lags.init[["tmax"   ]] <- data.frame(array(tmax.init   , dim=c(1, ens.hr)))
    lags.init[["tmin"   ]] <- data.frame(array(tmin.init   , dim=c(1, ens.hr)))
    lags.init[["precipf"]] <- data.frame(array(precipf.init, dim=c(1, ens.hr)))
    lags.init[["swdown" ]] <- data.frame(array(swdown.init , dim=c(1, ens.hr)))
    lags.init[["lwdown" ]] <- data.frame(array(lwdown.init , dim=c(1, ens.hr)))
    lags.init[["press"  ]] <- data.frame(array(press.init  , dim=c(1, ens.hr)))
    lags.init[["qair"   ]] <- data.frame(array(qair.init   , dim=c(1, ens.hr)))
    lags.init[["wind"   ]] <- data.frame(array(wind.init   , dim=c(1, ens.hr)))
  
  #for(y in years.sim){ # will need to add this for loop back in this maf
    dat.ens <- list() # a new list for each ensemble member as a new layer
    
    ens.sims <- list() # this will propogate that spread through each year, so instead of 
    # restarting every January 1, it will propogate those lag values
    # Create a list layer for each ensemble member
    nc.now <- ncdf4::nc_open(path.gcm)
    dat.yr <- data.frame(time    = ncdf4::ncvar_get(nc.now, "time"   ),
                         tmax    = ncdf4::ncvar_get(nc.now, "air_temperature_max"   ),
                         tmin    = ncdf4::ncvar_get(nc.now, "air_temperature_min"   ),
                         precipf = ncdf4::ncvar_get(nc.now, "precipitation_flux"),
                         swdown  = ncdf4::ncvar_get(nc.now, "surface_downwelling_shortwave_flux_in_air" ),
                         lwdown  = ncdf4::ncvar_get(nc.now, "surface_downwelling_longwave_flux_in_air" ),
                         press   = ncdf4::ncvar_get(nc.now, "air_pressure" ),
                         qair    = ncdf4::ncvar_get(nc.now, "specific_humidity"   ),
                         ewind    = ncdf4::ncvar_get(nc.now, "eastward_wind"   ),
                         nwind    = ncdf4::ncvar_get(nc.now, "northward_wind"   )
                         
    )
    dat.yr$wind = sqrt(dat.yr$ewind^2+dat.yr$nwind^2)
    ncdf4::nc_close(nc.now)
    
    # Do some stuff to get the right time variables for dat.yr
    dat.yr$year <- y
    dat.yr$date <- as.Date(dat.yr$time/86400, origin= paste0(y-1,"-12-31")) #doy needs to start at 1 here, not 0
    dat.yr$doy <- lubridate::yday(dat.yr$date)
    
    # Create the data frame for the "next" values
    dat.nxt <- dat.yr
    # Shift everyting up by a day to get the preview of the next day to get processed
    # Note: Because we work backwards through time, Jan 1 is the "next" day to get processed with Jan 2
    dat.nxt[2:(nrow(dat.nxt)), c("tmax", "tmin", "precipf", 
                                 "swdown", "lwdown", "press", "qair", "wind")] <- dat.nxt[1:(nrow(dat.nxt)-1),
                                                                                  c("tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
    # next day associates january 2nd with january 1
    #might need to flip this maf to +1 in dat.nxt[1:(nrow(dat.nxt))] might need to dat.nxt[2:(nrow(dat.nxt))]
    
    # Need to add in the "next" value for january (dec 31 of next year) 
    # Note: if we're past the end of our daily data, the best we can do is leave things as is (copy the last day's value)
    if(y>min(years)){
      
      nc.nxt <- ncdf4::nc_open(path.gcm)
      # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
      nxt.time <- ncdf4::ncvar_get(nc.nxt, "time")
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmax"   ] <- ncdf4::ncvar_get(nc.nxt, "air_temperature_max")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmin"   ] <- ncdf4::ncvar_get(nc.nxt, "air_temperature_min")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"precipf"] <- ncdf4::ncvar_get(nc.nxt, "precipitation_flux")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"swdown" ] <- ncdf4::ncvar_get(nc.nxt, "surface_downwelling_shortwave_flux_in_air")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"lwdown" ] <- ncdf4::ncvar_get(nc.nxt, "surface_downwelling_longwave_flux_in_air")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"press"  ] <- ncdf4::ncvar_get(nc.nxt, "air_pressure")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"qair"   ] <- ncdf4::ncvar_get(nc.nxt, "specific_humidity")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"ewind"  ] <- ncdf4::ncvar_get(nc.nxt, "eastward_wind")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"nwind"  ] <- ncdf4::ncvar_get(nc.nxt, "northward_wind")[length(nxt.time)]
      dat.nxt[dat.nxt$time==min(dat.nxt$time),"wind"   ] <- sqrt(dat.nxt$ewind^2 + dat.nxt$nwind^2)[length(nxt.time)]
      
      ncdf4::nc_close(nc.nxt)
    } 
    
    dat.ens <- data.frame(
                          year         =dat.yr $year   ,
                          doy          =dat.yr $doy    ,
                          date         =dat.yr $date   ,
                          tmax.day     =dat.yr $tmax   ,
                          tmin.day     =dat.yr $tmin   ,
                          precipf.day  =dat.yr $precipf,
                          swdown.day   =dat.yr $swdown ,
                          lwdown.day   =dat.yr $lwdown ,
                          press.day    =dat.yr $press  ,
                          qair.day     =dat.yr $qair   ,
                          wind.day     =dat.yr $wind   ,
                          next.tmax    =dat.nxt$tmax   ,
                          next.tmin    =dat.nxt$tmin   ,
                          next.precipf =dat.nxt$precipf,
                          next.swdown  =dat.nxt$swdown ,
                          next.lwdown  =dat.nxt$lwdown ,
                          next.press   =dat.nxt$press  ,
                          next.qair    =dat.nxt$qair   ,
                          next.wind    =dat.nxt$wind
    )
    
    #so start point should be start piont -1, so we changed this to 2018-12-31
    prev_year = y -1
    dat.ens$time.day <- as.numeric(difftime(dat.ens$date , paste0(y-1, "-12-31"), tz="GMT", units="day"))
    dat.ens <- merge(dat.ens, df.hour, all=T)
    
    dat.ens$date <- strptime(paste(dat.ens$year, dat.ens$doy, dat.ens$hour, sep="-"), "%Y-%j-%H", tz="GMT")
    dat.ens$time.hr <- as.numeric(difftime(dat.ens$date, paste0(y-1,"-12-31"), tz="GMT", units="hour")) #+ minute(dat.train$date)/60
    dat.ens <- dat.ens[order(dat.ens$time.hr),]
    
    # -----------------------------------
    # 2. Predict met vars for each ensemble member
    # Note: Using a loop for each ensemble member for now, but this will get 
    #       parallelized to speed it up soon, but we'll prototype in parallel
    # -----------------------------------
    
    ens.sims <- predict.subdaily(dat.ens, n.ens=ens.hr, path.model=file.path(lm.models.base), lags.list=NULL, lags.init=lags.init, dat.train=dat.train)

  # Set up the time dimension for this year
    hrs.now <- as.numeric(difftime(dat.ens$date, paste0(y,"-01-01"), tz="GMT", units="hour"))

    for(v in names(ens.sims)){
      lags.init[[v]] <- data.frame(ens.sims[[v]][length(ens.sims[[v]]),])
    }
      
      # -----------------------------------
      # Write each year for each ensemble member into its own .nc file
      # -----------------------------------
    ntime = nrow(dat.ens)
    lat <- ncdf4::ncdim_def(name = "latitude", units = "degree_north", vals = lat.in, create_dimvar = TRUE)
    lon <- ncdf4::ncdim_def(name = "longitude", units = "degree_east", vals = lon.in, create_dimvar = TRUE)
    time <- ncdf4::ncdim_def(name = "time", units = "sec", vals = (1:ntime) * 3600, 
                             create_dimvar = TRUE, unlim = TRUE)
    dim <- list(lat, lon, time)
    
    var.list <- list()
    for (j in seq_along(vars.info$name)){
      var.list[[j]] <- ncdf4::ncvar_def(name = as.character(vars.info$name.cf[j]), 
                                      units = as.character(vars.info$units[j]), 
                                      dim = dim, 
                                      missval = -9999, 
                                      verbose = verbose)
    }
      
    for (i in seq_len(ens.hr)){
      df <- data.frame(matrix(ncol =  length(vars.info$name), nrow = ntime))
      colnames(df) <- vars.info$name
      for (j in vars.info$name){
        dat.sim[[j]][["X1"]]
        e = paste0("X",i)
        df[[j]] = dat.sim[[j]][[e]]
      }
      
      df <- df[,c("tair", "precipf","swdown","lwdown","press","qair","wind")]
      colnames(df) = vars.info$name.cf
  
  
      rows <- 1
      dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
      
      loc.file <- file.path(outfolder, paste0(in.prefix, "_ens", 
                                              i, "_", y, ".nc"))
  
      loc <- ncdf4::nc_create(filename = loc.file, vars = var.list, verbose = verbose)
      
      for (j in vars.info$name.cf) {
        ncdf4::ncvar_put(nc = loc, varid = as.character(j), vals = df[[j]][seq_len(nrow(df))])
      }
      ncdf4::nc_close(loc)
    }
  }
}
      
    
    
    