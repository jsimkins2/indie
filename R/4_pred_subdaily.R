##' Create statistical models to predict subdaily meteorology from daily means
##' 
##' @title gen_subdaily_models
##' @param outfolder - directory where output file will be stored
##' @param in.path - path to model dataset you wish to temporally downscale
##' @param in.prefix - prefix of model dataset, i.e. if file is GFDL.CM3.rcp45.r1i1p1.2006 the prefix is "GFDL.CM3.rcp45.r1i1p1"
##' @param dat.train_file - location of train_data file
##' @param start_date
##' @param end_date
##' @param tdf_file - temporal_downscale_functions.R filepath that can be sourced i.e. "~/scripts/temporal_downscale_functions.R"
##' @param td_file - temporal_downscale.R filepath that can be sourced i.e. "~/scripts/temporal_downscale.R"
##' @param site.name - name of location i.e. Ameriflux Site as character string
##' @param site.lat - latitude as numeric
##' @param site.lon - longitude as numeric
##' @param ens.hr - integer selecting number of hourly ensemble members 
##' @param n.day - 1 # Number of daily ensemble members to process, must be 1 because not passing in daily ensemble
##' @param years.sim - NULL
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

pred_subdaily_met <- function(outfolder, in.path, in.prefix, dat.train_file, start_date, end_date, 
                              tdf_file, td_file, site.name, site.lat, site.lon,
                              years.sim = NULL, cores.max = 12,  
                              ens.hr = 3, n.day = 1, resids = FALSE, 
                              parallel = FALSE, n.cores=NULL, day.window) {
  
  outfolder = "~"
  in.path = "~_site_0-2"
  in.prefix = "MACA.IPSL-CM5A-LR.rcp85.r1i1p1"
  dat.train_file = "US-WCr_train_data"
  start_date = "2019-01-01"
  end_date = "2022-12-31"
  tdf_file = "Temp_downscale/Scripts/temp_dwnsc_functions.R"
  td_file = "Temp_downscale/Scripts/temp_dwnsc.R"
  site.name = "US-WCr"
  site.lat = 45
  site.lon = -90
  day.window = 5
  years.sim = NULL
  cores.max = 12
  ens.hr = 3
  n.day = 1
  resids = FALSE
  parallel = FALSE
  n.cores = NULL
  lm.models.base = "sf_scratch/US-WCr/"
  
  years = seq(lubridate::year(start_date), lubridate::year(end_date))

  library(ncdf4)
  library(mgcv)
  library(MASS)
  library(lubridate)
  library(ggplot2)
  library(stringr)
  library(tictoc)
  library(parallel)

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
  # Make a few dimensions we can use
  dimY <- ncdim_def( "lon", units="degrees", longname="latitude", vals=site.lat )
  dimX <- ncdim_def( "lat", units="degrees", longname="longitude", vals=site.lon )
  # -----------------------------------
  
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
  
      
      nc.now <- nc_open(path.gcm)
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
      
      ens.sims <- list() # this will propogate that spread through each year, so instead of 
      # restarting every January 1, it will propogate those lag values
      # Create a list layer for each ensemble member
        nc.now <- nc_open(path.gcm)
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
        dat.yr$lwdown = 250
        dat.yr$press = 10000
        nc_close(nc.now)
        
        # Do some stuff to get the right time variables for dat.yr
        dat.yr$year <- y
        dat.yr$date <- as.Date(dat.yr$time/86400, origin= paste0(y,"-12-31")) #doy needs to start at 1 here, not 0
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
        if(y>min(years)){
          
          #file.next needs to be +1 here so we go forward in years
          
          nc.nxt <- nc_open(path.gcm)
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
          dat.nxt$lwdown = 250
          dat.nxt$press = 10000
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
        dat.ens$time.day <- as.numeric(difftime(dat.ens$date, paste0(y, "-12-31"), tz="GMT", units="day"))
        dat.ens <- merge(dat.ens, df.hour, all=T)
        
        # A very ugly hack way to get the minutes in there, yeah don't really need this maf,
        #can probably get rid of this and anywhere else that gets a minute aiight
        #dat.ens$minute <- abs(dat.ens$hour-(round(dat.ens$hour, 0)))*60
        
        dat.ens$date <- strptime(paste(dat.ens$year, dat.ens$doy, dat.ens$hour, sep="-"), "%Y-%j-%H", tz="GMT")
        dat.ens$time.hr <- as.numeric(difftime(dat.ens$date, paste0(y,"-12-31"), tz="GMT", units="hour")) #+ minute(dat.train$date)/60
        dat.ens <- dat.ens[order(dat.ens$time.hr),]
        
        ens.sims <- predict.subdaily(dat.ens, n.ens=ens.hr, path.model=file.path(lm.models.base), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
        # End ensembles setup
      
      # Set up the time dimension for this year
      hrs.now <- as.numeric(difftime(dat.ens$date, paste0(y,"-01-01"), tz="GMT", units="hour"))
      dim.t <- ncdim_def(name = "time",
                         units = paste0("hours since ", y, "-01-01 00:00:00:"),
                         vals = hrs.now, # calculating the number of months in this run
                         calendar = "standard", unlim = TRUE)
      
      
      # -----------------------------------
      # 2. Predict met vars for each ensemble member
      # Note: Using a loop for each ensemble member for now, but this will get 
      #       parallelized to speed it up soon, but we'll prototype in parallel
      # -----------------------------------
  
        for(v in names(ens.sims)){
          lags.init[[v]] <- data.frame(ens.sims[[v]][length(ens.sims[[v]]),])
        }
        
        # -----------------------------------
        # Write each year for each ensemble member into its own .nc file
        # -----------------------------------
        for(i in 1:ens.hr){
          ens.name <- paste0(site.name, "_", in.prefix, "_hr", "_", str_pad(e, 3, pad=0), "-", str_pad(i, 3, pad=0))
          
          if(!dir.exists(file.path(path.out, ens.name))) dir.create(file.path(path.out, ens.name), recursive=T)
          # path.out <- file.path(dat.base, GCM, "1hr")
          
          var.list <- list()
          dat.list <- list()
          for(v in names(ens.sims[[paste0("X", e)]])){
            var.cf = vars.info[vars.info$name==v, "name.cf"]
            var.list[[v]] <- ncvar_def(v, units=paste(vars.info[vars.info$name==v, "units"]), dim=list(dimX, dimY, dim.t), longname=paste(vars.info[vars.info$name==v, "longname"]))
            dat.list[[v]] <- array(ens.sims[[paste0("X", e)]][[v]][,i], dim=c(1,1,length(hrs.now)))
          }
          
          # Naming convention: [SITE]_[GCM]_1hr_[bias_ens_member]-[subday_ens_member]_[YEAR].nc
          nc <- nc_create(file.path(path.out, ens.name, paste0(ens.name,"_", str_pad(y, 4, pad=0), ".nc")), var.list)
          for(v in 1:length(var.list)) {
            ncvar_put(nc, var.list[[v]], dat.list[[v]])
          }
          nc_close(nc)    
        }
        # -----------------------------------
        
      } # End ensemble member prediction for 1 year
      # -----------------------------------
    
      # need to figure out a way to save these by ensemble member, and have each ensemble member work in conjunction
      reso_len = 8760
      y = year
      
      for (i in seq_len(ens.hr)){
        df <- data.frame(matrix(ncol =  length(vars.info$name), nrow = reso_len))
        colnames(df) <- vars.info$name
        for (j in vars.info$name){
          dat.sim[[j]][["X1"]]
          e = paste0("X",i)
          df[[j]] = dat.sim[[j]][[e]]
        }
        assign(paste0("ens",i),df)
      }
        
      
      rows <- 1
      dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
      
      loc.file <- file.path(outfolder, paste0(site.name,in.prefix, "ens", 
                                              i, "_", y, ".nc"))
      
      loc <- ncdf4::nc_create(filename = loc.file, vars = , verbose = verbose)
      for (j in seq_along(var$CF.name)) {
        ncdf4::ncvar_put(nc = loc, varid = as.character(var$CF.name[j]), vals = downscaled.met[[j]])
      }
      ncdf4::nc_close(loc)
      
      results[[e]] <- data.frame(file = loc.file, 
                                 host = rep(PEcAn.utils::fqdn(),rows), 
                                 mimetype = rep("application/x-netcdf",rows), 
                                 formatname = rep("CF Meteorology",rows),
                                 startdate = paste0(year, "-01-01 00:00:00", tz = "UTC"), 
                                 enddate = paste0(year, "-12-31 23:59:59", tz = "UTC"),
                                 dbfile.name = paste0(source_name, ".dwnsc.ens"),
                                 stringsAsFactors = FALSE)
      
      
    
    return(invisible(results))    
}
      
    
    
    