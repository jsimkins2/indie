##' Predict subdaily meteorology from coarse datasets from training dataset linear regression models
##' 
##' @param outfolder - directory where output will be stored
##' @param in.path - path of coarse model (e.g. GCM output)
##' @param in.prefix - prefix of model string as character (e.g. IPSL.r1i1p1.rcp85)
##' @param start_date 
##' @param end_date
##' @param ens.hr - # Number of hourly ensemble members to create
##' @param dat.train - training dataset created by CF2traindata.R
##' @param tdf_filepath - temporal_downscale_functions.R path
##' @param tdp_filepath - tem_dwnsc_pred.R path 
##' @param n.beta - number of betas to save from linear regression model
##' @param resids - logical stating whether to pass on residual data or not
##' @param parallel - logical stating whether to run temporal_downscale_functions.R in parallel 
##' @param n.cores - deals with parallelization
##' @param day.window - integer specifying number of days around a particular day you want to use data from for that 
##'                     specific hours coefficients
##' @param overwrite
##' @param verbose

##' @author Christy Rollinson, James Simkins

pred_subdaily_met <- function(outfolder, in.path, in.prefix, start_date, end_date,
                              ens.hr=3, dat.train, tdf_filepath, tdp_filepath, n.beta,
                              resids=F, parallel=F, n.cores=NULL, day.window,
                              overwrite = TRUE, verbose = FALSE){
  
  
  start_date <- as.POSIXlt(start_date, tz = "UTC")
  end_date <- as.POSIXlt(end_date, tz = "UTC")
  start_year <- lubridate::year(start_date)
  end_year   <- lubridate::year(end_date)
  year_seq <- seq(start_year,end_year,1)
  
  library(ncdf4)
  library(mgcv)
  library(MASS)
  library(lubridate)
  library(ggplot2)
  library(stringr)
  library(tictoc)
  library(parallel)
  # library(tictoc)
   
  
  # Load the scripts that do all the heavy lifting
  source(tdf_filepath)
  source(tdp_filepath)

  dat.train <- read.csv(dat.train)
  
  # Hard-coding numbers for PARKFALLS
  site.name <- dat.train$site.name
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
  # Make a few dimensions we can use
  dimY <- ncdim_def( "lon", units="degrees", longname="latitude", vals=site.lat )
  dimX <- ncdim_def( "lat", units="degrees", longname="longitude", vals=site.lon )
  # -----------------------------------

  # Set & create the output directory
  path.out <- outfolder
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # -----------------------------------
  # 1. Format output so all ensemble members can be run at once
  # NOTE: Need to start with the last and work to the first
  # -----------------------------------
  # Figure out which daily ensembles we should pull from
  day.dirs <- dir(path.gcm,  paste0(site.name, "_", GCM, "_day_")) #directory of input data aka maca
  ens.day.all <- substr(day.dirs, nchar(day.dirs)-2, nchar(day.dirs))
  
  hrs.dir <- dir(file.path(path.out),  paste0(site.name, "_", GCM, "_1hr_")) #this is the output path
  ens.day.done <- substr(hrs.dir, nchar(hrs.dir)-6, nchar(hrs.dir)-4)
  
  ens.list <- ens.day.all[!(ens.day.all %in% ens.day.done)] #list of files so aka multiple years
  
  seed <- seed.vec[length(hrs.dir)+1] # This makes sure that if we add ensemble members, it gets a new, but reproducible seed
  set.seed(seed)
  
  
  
  # Initialize the lags
  for (y in year_seq){
    lags.init <- list() # Need to initialize lags first outside of the loop and then they'll get updated internally
      
      nc.now <- ncdf4::nc_open(paste0(path.in,prefix.in,year_seq[y]))
      # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
      nc.time <- ncvar_get(nc.now, "time")
      tmax.init    <- ncvar_get(nc.now, "air_temperature_max")[length(nc.time)]
      tmin.init    <- ncvar_get(nc.now, "air_temperature_min")[length(nc.time)]
      precipf.init <- ncvar_get(nc.now, "precipitation_flux")[length(nc.time)]
      swdown.init  <- ncvar_get(nc.now, "surface_downwelling_shortwave_flux_in_air")[length(nc.time)]
      lwdown.init  <- ncvar_get(nc.now, "surface_downwelling_longwave_flux_in_air")[length(nc.time)]
      press.init   <- ncvar_get(nc.now, "air_pressure")[length(nc.time)]
      qair.init    <- ncvar_get(nc.now, "specific_humidity")[length(nc.time)]
      nwind.init    <- ncvar_get(nc.now, "northward_wind")[length(nc.time)]
      ewind.init <- ncvar_get(nc.now, "eastward_wind")[length(nc.time)]
      wind.init <- sqrt(nwind.init^2 + ewind.init^2)
      nc_close(nc.now)
      
      lwdown.init = 250
      press.init = 10000
      
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
      
      #ens.sims <- list() # this will propogate that spread through each year, so instead of 
      # restarting every January 1, it will propogate those lag values
      
      # Create a list layer for each ensemble member
      # NOTE: NEED TO CHECK THE TIME STAMPS HERE TO MAKE SURE WE'RE PULLING AND ADDING IN THE RIGHT PLACE
        nc.now <- ncdf4::nc_open(paste0(path.in,prefix.in,year_seq[y]))
        dat.yr <- data.frame(time    = ncvar_get(nc.now, "time"   ),
                             tmax    = ncvar_get(nc.now, "air_temperature_max"   ),
                             tmin    = ncvar_get(nc.now, "air_temperature_min"   ),
                             precipf = ncvar_get(nc.now, "precipitation_flux"),
                             swdown  = ncvar_get(nc.now, "surface_downwelling_shortwave_flux_in_air" ),
                             lwdown  = ncvar_get(nc.now, "surface_downwelling_longwave_flux_in_air" ),
                             press   = ncvar_get(nc.now, "air_pressure" ),
                             qair    = ncvar_get(nc.now, "specific_humidity"   ),
                             ewind    = ncvar_get(nc.now, "eastward_wind"   ),
                             nwind    = ncvar_get(nc.now, "northward_wind"   )
                             
        )
        dat.yr$wind = sqrt(dat.yr$ewind^2+dat.yr$nwind^2)
        #dat.yr$lwdown = 250 
        #dat.yr$press = 10000
        nc_close(nc.now)
        
        # Do some stuff to get the right time variables for dat.yr
        dat.yr$year <- y
        dat.yr$date <- as.Date(dat.yr$time/86400, origin=paste0(y,"-12-31")) #doy needs to start at 1 here, not 0
        dat.yr$doy <- yday(dat.yr$date)
        
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
        if(y>min(years.sim)){
          
          #file.next needs to be +1 here so we go forward in years
          
          nc.nxt <- ncdf4::nc_open(paste0(path.in,prefix.in,year_seq[y]))
          # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
          nxt.time <- ncvar_get(nc.nxt, "time")
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmax"   ] <- ncvar_get(nc.nxt, "air_temperature_max")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmin"   ] <- ncvar_get(nc.nxt, "air_temperature_min")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"precipf"] <- ncvar_get(nc.nxt, "precipitation_flux")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"swdown" ] <- ncvar_get(nc.nxt, "surface_downwelling_shortwave_flux_in_air")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"lwdown" ] <- ncvar_get(nc.nxt, "surface_downwelling_longwave_flux_in_air")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"press"  ] <- ncvar_get(nc.nxt, "air_pressure")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"qair"   ] <- ncvar_get(nc.nxt, "specific_humidity")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"ewind"   ] <- ncvar_get(nc.nxt, "eastward_wind")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"nwind"   ] <- ncvar_get(nc.nxt, "northward_wind")[length(nxt.time)]
          dat.nxt[dat.nxt$time==min(dat.nxt$time),"wind"   ] <- sqrt(dat.nxt$ewind^2 + dat.nxt$nwind^2)[length(nxt.time)]
          
          nc_close(nc.nxt)
          #dat.nxt$lwdown = 250
          #dat.nxt$press = 10000
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
        dat.ens$time.day <- as.numeric(difftime(dat.ens$date, paste0(y,"-12-31"), tz="GMT", units="day"))
        dat.ens <- merge(dat.ens, df.hour, all=T)
        
        # A very ugly hack way to get the minutes in there, yeah don't really need this maf,
        #can probably get rid of this and anywhere else that gets a minute aiight
        #dat.ens$minute <- abs(dat.ens$hour-(round(dat.ens$hour, 0)))*60
        
        # timestep <- max(dat.ens[[paste0("X", e)]]$minute)/60 # Calculate the timestep to make the date function work
        dat.ens$date <- strptime(paste(dat.ens$year, dat.ens$doy, dat.ens$hour, sep="-"), "%Y-%j-%H", tz="GMT")
        dat.ens$time.hr <- as.numeric(difftime(dat.ens$date, paste0(y,"-12-31"), tz="GMT", units="hour")) #+ minute(dat.train$date)/60
        dat.ens <- dat.ens[order(dat.ens$time.hr),]
        
        ens.sims <- predict.subdaily(dat.ens, n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
        # End ensembles setup
      
      # Set up the time dimension for this year
      hrs.now <- as.numeric(difftime(dat.ens$date, paste0(y,"-01-01"), tz="GMT", units="hour"))
      dim.t <- ncdim_def(name = "time",
                         units = paste0("hours since ",y,"-01-01 00:00:00:"),
                         vals = hrs.now, # calculating the number of months in this run
                         calendar = "standard", unlim = TRUE)
      
      
      # -----------------------------------
      # 2. Predict met vars for each ensemble member
      # Note: Using a loop for each ensemble member for now, but this will get 
      #       parallelized to speed it up soon, but we'll prototype in parallel
      # -----------------------------------
      #cores.use <- min(cores.max, length(dat.ens))
      #ens.sims  <- mclapply(dat.ens, predict.subdaily, mc.cores=cores.use, n.ens=ens.hr, path.model="Ameri_downscale/", lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
      
      # ens.sims <- predict.subdaily(dat.ens[[1]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.init=lags.init[[1]], dat.train=dat.train)
      
    
      
    
        # Update the initial lags for next year
        for(v in names(ens.sims)){
          lags.init[[v]] <- data.frame(ens.sims[[v]][length(ens.sims[[v]]),])
        }
    
        # -----------------------------------
      # -----------------------------------
    
      # need to figure out a way to save these by ensemble member, and have each ensemble member work in conjunction
      reso_len = 8760
      
      for (i in seq_len(ens.hr)){
        df <- data.frame(matrix(ncol =  length(vars.info$name), nrow = reso_len))
        colnames(df) <- vars.info$name
        for (j in vars.info$name){
          dat.sim[["tair"]][["X1"]]
          e = paste0("X",i)
          df[[j]] = dat.sim[[j]][[e]]
        }
        assign(paste0("ens",i),df)
      }
        
      
      rows <- 1
      dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
      
      loc.file <- file.path(outfolder, paste0(site.name,GCM, "ens", 
                                              i, "_", y, ".nc"))
      
      loc <- ncdf4::nc_create(filename = loc.file, vars = train.list, verbose = verbose)
      for (j in seq_along(var$CF.name)) {
        ncdf4::ncvar_put(nc = loc, varid = as.character(var$CF.name[j]), vals = downscaled.met[[j]])
      }
      ncdf4::nc_close(loc)
      
      results[[e]] <- data.frame(file = loc.file, 
                                 host = rep(PEcAn.utils::fqdn(),rows), 
                                 mimetype = rep("application/x-netcdf",rows), 
                                 formatname = rep("CF Meteorology",rows),
                                 startdate = paste0(y, "-01-01 00:00:00", tz = "UTC"), 
                                 enddate = paste0(y, "-12-31 23:59:59", tz = "UTC"),
                                 dbfile.name = paste0(source_name, ".dwnsc.ens"),
                                 stringsAsFactors = FALSE)
      
      }
    
    return(invisible(results))    
    }
  }
}
    
    
    
    