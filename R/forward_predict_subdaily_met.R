# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Create statistical models to predict subdaily meteorology from daily means
# Creator: Christy Rollinson, 28 October 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Apply the statistical models from step 3 to convert the daily, bias-corrected met 
# files from step 2 (daily means) and predict subdaily values.  This gets done by 
# filtering backwards in time starting with the present (where the trianing data is).
#
# There are ways to improve this and speed it up, but hopefully this works for now.
# We whould also probably think about applying this filter approach to the bias-
# correction step to avoid abrupt and unreasonable jumps in climate.
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Load libraries, set up file paths, etc
# ----- Loop through by ensemble member ----------
#    1. Load and format prediction data (1 file from step 2)
#       1.1 Load output file from bias correction (bias ensemble member)
#       1.2 select year we're working with
#    ----- Loop through by year ----------
#      2. Predict subdaily values for whole year, filtering backwards in time
#      3. Write annual output into .nc files 
#         - separate file for each year/ensemle member; 
#         - all met vars in one annual file (similar to pecan met structure)
#    ----- recycle steps 2 & 3 for all years in file ----------
# ----- recycle step 1 for all files for ensemble member ----------
# -----------------------------------


# -----------------------------------
# 0. Load libraries, set up file paths, etc
# -----------------------------------
# Script to prototype temporal downscaling
library(ncdf4)
library(mgcv)
library(MASS)
library(lubridate)
library(ggplot2)
library(stringr)
library(tictoc)
library(parallel)
# library(tictoc)
rm(list=ls())

wd.base <- "~/Christy_code"
# wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble/"
setwd(wd.base)

# Load the scripts that do all the heavy lifting
source("temporal_downscale.R")
source("temporal_downscale_functions.R")


dat.base <- "~/Christy_code/Christy_code/maca_data_site_0-1"
# dat.base <- "~/Desktop/met_ensembles/PARKFALLS/"
# dat.base <- "~/Desktop/met_bias_day/data/met_ensembles/PARKFALLS/"

dat.train <- read.csv(file.path(wd.base, "PFa_training_data_com"))

# Hard-coding numbers for PARKFALLS
site.name="PARKFALLS"
site.lat=45.95
site.lon=-90.27

GCM = c("IPSL")
# GCM.list = "MIROC-ESM"
ens.hr  <- 3 # Number of hourly ensemble members to create
n.day <- 10 # Number of daily ensemble members to process
yrs.plot <- c(2020)
# years.sim=2015:1900
years.sim=NULL
cores.max = 12

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

# NOTE: all precip needs to be converted precip back to kg/m2/s from kg/m2/day
# This gets done when formatting things for downscaling

  # Set the directory where the output is & load the file
  path.gcm <- dat.base

  # Set & create the output directory
  path.out <- file.path(dat.base, GCM, "1hr")
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # -----------------------------------
  # 1. Format output so all ensemble members can be run at once
  # NOTE: Need to start with the last and work to the first
  # -----------------------------------

  seed <- seed.vec[length(hrs.dir)+1] # This makes sure that if we add ensemble members, it gets a new, but reproducible seed
  set.seed(seed)

  
  # Initialize the lags
  lags.init <- list() # Need to initialize lags first outside of the loop and then they'll get updated internally

    GCM_file = file.path(path.gcm, "MACA.IPSL-CM5A-LR.rcp85.r1i1p1.2020.nc")
    nc.now <- nc_open(GCM_file)
    # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
    nc.time <- ncvar_get(nc.now, "time")
    tmax.init    <- ncvar_get(nc.now, "air_temperature_max")[length(nc.time)]
    tmin.init    <- ncvar_get(nc.now, "air_temperature_min")[length(nc.time)]
    precipf.init <- ncvar_get(nc.now, "precipitation_flux")[length(nc.time)]
    swdown.init  <- ncvar_get(nc.now, "surface_downwelling_shortwave_flux_in_air")[length(nc.time)]
    lwdown.init  <- ncvar_get(nc.now, "surface_downwelling_longwave_flux_in_air")[length(nc.time)]
    press.init   <- ncvar_get(nc.now, "air_pressure")[length(nc.time)]
    qair.init    <- ncvar_get(nc.now, "specific_humidity")[length(nc.time)]
    ewind.init   <- ncvar_get(nc.now, "eastward_wind")[length(nc.time)]
    nwind.init   <- ncvar_get(nc.now, "northward_wind")[length(nc.time)]
    wind.init <- sqrt(ewind.init^2+nwind.init^2)
    nc_close(nc.now)
    
    lags.init[[paste0("X")]][["tair"   ]] <- data.frame(array(mean(c(tmax.init, tmin.init)), dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["tmax"   ]] <- data.frame(array(tmax.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["tmin"   ]] <- data.frame(array(tmin.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["precipf"]] <- data.frame(array(precipf.init, dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["swdown" ]] <- data.frame(array(swdown.init , dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["lwdown" ]] <- data.frame(array(lwdown.init , dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["press"  ]] <- data.frame(array(press.init  , dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["qair"   ]] <- data.frame(array(qair.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X")]][["wind"   ]] <- data.frame(array(wind.init   , dim=c(1, ens.hr)))
  
  years.sim = 2019
  y = 2019
    dat.ens <- list() # a new list for each ensemble member as a new layer
    df.hour <- data.frame(hour=0:23)
    
    # Create a list layer for each ensemble member
    # NOTE: NEED TO CHECK THE TIME STAMPS HERE TO MAKE SURE WE'RE PULLING AND ADDING IN THE RIGHT PLACE

      nc.now <- nc_open(GCM_file)
      nwind    = ncvar_get(nc.now, "northward_wind"   )
      ewind   = ncvar_get(nc.now, "eastward_wind"      )
      dat.yr <- data.frame(time    = ncvar_get(nc.now, "time"   ),
                           tmax    = ncvar_get(nc.now, "air_temperature_max"   ),
                           tmin    = ncvar_get(nc.now, "air_temperature_min"   ),
                           precipf = ncvar_get(nc.now, "precipitation_flux"),
                           swdown  = ncvar_get(nc.now, "surface_downwelling_shortwave_flux_in_air" ),
                           lwdown  = ncvar_get(nc.now, "surface_downwelling_longwave_flux_in_air" ),
                           press   = ncvar_get(nc.now, "air_pressure"  ),
                           qair    = ncvar_get(nc.now, "specific_humidity"   ),
                           wind = sqrt(nwind^2+ewind^2))
    
      nc_close(nc.now)
      
      # Do some stuff to get the right time variables for dat.yr
      dat.yr$year <- y
      dat.yr$date <- as.Date(dat.yr$time, origin="2019-01-01")
      dat.yr$doy <- yday(dat.yr$date)
      
      # Create the data frame for the "next" values
      dat.nxt <- dat.yr
      # Shift everyting up by a day to get the preview of the next day to get processed
      # Note: Because we work backwards through time, Jan 1 is the "next" day to get processed with Jan 2
      dat.nxt[2:(nrow(dat.nxt)), c("tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind")] <- dat.nxt[1:(nrow(dat.nxt)-1), c("tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
      
      # Need to add in the "next" value for january (dec 31 of next year) 
      # Note: if we're past the end of our daily data, the best we can do is leave things as is (copy the last day's value)
        
        nc.nxt <- nc_open(GCM_file)
        # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
        nxt.time <- ncvar_get(nc.nxt, "time")
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmax"   ] <- ncvar_get(nc.nxt, "air_temperature_max")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmin"   ] <- ncvar_get(nc.nxt, "air_temperature_min")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"precipf"] <- ncvar_get(nc.nxt, "precipitation_flux")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"swdown" ] <- ncvar_get(nc.nxt, "surface_downwelling_shortwave_flux_in_air")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"lwdown" ] <- ncvar_get(nc.nxt, "surface_downwelling_longwave_flux_in_air")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"press"  ] <- ncvar_get(nc.nxt, "air_pressure")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"qair"   ] <- ncvar_get(nc.nxt, "specific_humidity")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"nwind"  ] <- ncvar_get(nc.nxt, "northward_wind")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"ewind"  ] <- ncvar_get(nc.nxt, "northward_wind")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"wind"   ] <- sqrt(dat.nxt$nwind^2 + dat.nxt$ewind^2)[length(nxt.time)]
        
        nc_close(nc.nxt)
      
      dat.ens[[paste0("X")]] <- data.frame(ens.day      =as.factor(paste0("X")),
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
      
      dat.ens[[paste0("X")]]$time.day <- as.numeric(difftime(dat.ens[[paste0("X")]]$date, "2019-01-01", tz="GMT", units="day"))
      dat.ens[[paste0("X")]] <- merge(dat.ens[[paste0("X")]], df.hour, all=T)
      
      dat.ens[[paste0("X")]]$date <- strptime(paste(dat.ens[[paste0("X")]]$year, dat.ens[[paste0("X")]]$doy, dat.ens[[paste0("X")]]$hour, sep="-"), "%Y-%j-%H", tz="GMT")
      dat.ens[[paste0("X")]]$time.hr <- as.numeric(difftime(dat.ens[[paste0("X")]]$date, "2019-01-01", tz="GMT", units="hour"))
      dat.ens[[paste0("X")]] <- dat.ens[[paste0("X")]][order(dat.ens[[paste0("X")]]$time.hr),]
    
    # Set up the time dimension for this year
    hrs.now <- as.numeric(difftime(dat.ens[[paste0("X")]]$date, "2019-01-01", tz="GMT", units="hour"))
    dim.t <- ncdim_def(name = "time",
                       units = paste0("hours since 2019-01-01 00:00:00:"),
                       vals = hrs.now, # calculating the number of months in this run
                       calendar = "standard", unlim = TRUE)
    
    
    # -----------------------------------
    # 2. Predict met vars for each ensemble member
    # Note: Using a loop for each ensemble member for now, but this will get 
    #       parallelized to speed it up soon, but we'll prototype in parallel
    # -----------------------------------
    cores.use <- min(cores.max, length(dat.ens))
    model.base <- "~/Christy_code/Christy_code/"
    ens.sims  <- mclapply(dat.ens, predict.subdaily, mc.cores=cores.use, n.ens=ens.hr, path.model=file.path(model.base, "mod_out"), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
    
    # If this is one of our designated QAQC years, makes some graphs
    # Now doing this for the whole GCM
    if(y %in% yrs.plot){
      dat.plot <- data.frame()
      ens.plot <- list()
      for(i in names(ens.sims[[1]])){
        ens.plot[[i]] <- data.frame(matrix(nrow=nrow(ens.sims[[1]][[i]]), ncol=0))
      }
      for(e in names(ens.sims)){
        dat.plot <- rbind(dat.plot, dat.ens[[e]])
        for(i in names(ens.sims[[e]])){
          ens.plot[[i]] <- cbind(ens.plot[[i]], ens.sims[[e]][[i]])
        }
      }
      day.name <- paste0(site.name, "_", GCM, "_1hr")
      fig.ens <- file.path(path.out, "subdaily_qaqc", day.name)
      if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
      for(v in names(ens.plot)){
        graph.predict(dat.mod=dat.plot, dat.ens=ens.plot, var=v, fig.dir=fig.ens)
      }
    }
    
    
    # ens.sims <- list()
    for(e in unique(ens.day)){
      # # Do the prediction
      # ens.sims[[paste0("X", e)]] <- predict.subdaily(dat.mod=dat.ens[[paste0("X", e)]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
      # qair.max <- quantile(as.matrix(ens.sims[[paste0("X", e)]]$qair[,c(1:ens.hr)]), 0.99)
      # ens.sims[[paste0("X", e)]]$qair[ens.sims[[paste0("X", e)]]$qair>qair.max] <- qair.max
      
      # Update the initial lags for next year
      for(v in names(ens.sims[[paste0("X", e)]])){
        lags.init[[paste0("X",e)]][[v]] <- data.frame(ens.sims[[paste0("X", e)]][[v]][length(ens.sims[[paste0("X", e)]][[v]]),])
      }
      
      # -----------------------------------
      # Write each year for each ensemble member into its own .nc file
      # -----------------------------------
      for(i in 1:ens.hr){
        ens.name <- paste0(site.name, "_", GCM, "_1hr_", str_pad(e, 3, pad=0), "-", str_pad(i, 3, pad=0))
        
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
  } # End Year Loop
  # -----------------------------------
  
  # Do some clean-up to save space
  # dir.compress <- dir(path.out, GCM)
  
  # setwd(path.out)
  # for(ens in dir.compress){
  #   system(paste0("tar -jcvf ", ens, ".tar.bz2 ", ens)) # Compress the folder
  #   system(paste0("rm -rf ", ens)) # remove the uncompressed folder
  # }
  # setwd(wd.base)
  toc()
} # End GCM loop