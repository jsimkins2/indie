# Set progress bar
pb.index=1
pb <- txtProgressBar(min=1, max=6, style=3)

predict.subdaily <- function(dat.mod, n.ens, path.model, lags.list=NULL, lags.init=NULL, dat.train){
  
  # dat.mod    = data to be predicted at the time step of the training data
  # n.ens      = number of hourly ensemble members to generate
  # path.model = path to where the training model & betas is stored
  # lags.list  = a list with layers of same name dat.sim and n=n.ens that provide the initial lags; 
  #              used if entering the function from a parallel apply function
  # lags.init  = a data frame of initialization paramters to match the data in dat.mod
  # dat.train  = the training data used to fit the model; needed for night/day in surface_downwelling_shortwave_flux_in_air
  
  # --------------------------------
  # Models each variable separately at the moment rather than a generalized eqauiton that could potentially be parallilizeable
  # --------------------------------
  # #    2.1 air_temperature (air temperature)
  # #    2.2 precipitation.flux (preciptiation, water equivalent)
  # #    2.3 surface_downwelling_shortwave_flux_in_air (downwelling shortwave radiation)
  # #    2.4 surface_downwelling_longwave_flux_in_air (downwelling longwave radiation)
  # #    2.5 air_pressure (surface air_pressureure)
  # #    2.6 specific_humidity (specific humidity)
  # #    2.7 wind_speed (surface wind_speed speed)
  # --------------------------------# 
  # Load libraries
  #library(ncdf4)
  #library(mgcv)
  #library(MASS)
  #library(ggplot2)
  
  setTxtProgressBar(pb, pb.index)

  # Figure out if we need to extract the approrpiate 
  if(is.null(lags.init)){
    lags.init <- lags.list[[unique(dat.mod$ens.day)]]
  }
  
  # Set up the ensemble members in a list so the uncertainty can be propogated
  dat.sim <- list() 
  
  # DOY indexing is now off from the original; fix by subtracting 1
  # dat.mod$doy = dat.mod$doy-1
  # ------------------------------------------
  # Modeling surface_downwelling_shortwave_flux_in_air 
  # Note: this can be generalized to just run by DOY for all years at once since there's no memory in the system
  # ------------------------------------------
  {
    # Load the meta info for the betas
    betas.surface_downwelling_shortwave_flux_in_air <- ncdf4::nc_open(file.path(path.model, "surface_downwelling_shortwave_flux_in_air", "betas_surface_downwelling_shortwave_flux_in_air_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.surface_downwelling_shortwave_flux_in_air, "1"))
    ncdf4::nc_close(betas.surface_downwelling_shortwave_flux_in_air)
    
    dat.sim[["surface_downwelling_shortwave_flux_in_air"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      # For surface_downwelling_shortwave_flux_in_air, we only want to model daylight hours -- make sure this matches what's in the surface_downwelling_shortwave_flux_in_air function
      day.now = unique(dat.mod[dat.mod$time.day==i, "doy"])
      
      # Use the training data to figure out night/day
      hrs.day = unique(dat.train[dat.train$doy==day.now & dat.train$surface_downwelling_shortwave_flux_in_air>quantile(dat.train[dat.train$surface_downwelling_shortwave_flux_in_air>0,"surface_downwelling_shortwave_flux_in_air"], 0.05), "hour"])
      
      rows.now = which(dat.mod$time.day==i)
      rows.mod = which(dat.mod$time.day==i & dat.mod$hour %in% hrs.day)
      
      dat.temp <- dat.mod[rows.mod,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      
      dat.temp$surface_downwelling_shortwave_flux_in_air = 99999 # Dummy value so there's a column
      # day.now = unique(dat.temp$doy)
      
      # Load the saved model
      load(file.path(path.model, "surface_downwelling_shortwave_flux_in_air", paste0("model_surface_downwelling_shortwave_flux_in_air_", day.now, ".Rdata")))
      mod.surface_downwelling_shortwave_flux_in_air.doy <- mod.save
      rm(mod.save)
      
      # Pull coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.surface_downwelling_shortwave_flux_in_air <- ncdf4::nc_open(file.path(path.model, "surface_downwelling_shortwave_flux_in_air", paste0("betas_surface_downwelling_shortwave_flux_in_air_", day.now, ".nc")))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.surface_downwelling_shortwave_flux_in_air, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.surface_downwelling_shortwave_flux_in_air)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.surface_downwelling_shortwave_flux_in_air.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      dat.pred[dat.pred<0] <- 0
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$surface_downwelling_shortwave_flux_in_air), replace=T)
      for(j in 1:ncol(dat.sim$surface_downwelling_shortwave_flux_in_air)){
        dat.sim[["surface_downwelling_shortwave_flux_in_air"]][rows.mod,j] <- dat.pred[,cols.prop[j]]
      }
      
      # For night time hours, value shoudl be 0
      dat.sim[["surface_downwelling_shortwave_flux_in_air"]][rows.now[!rows.now %in% rows.mod],] <- 0
      rm(mod.surface_downwelling_shortwave_flux_in_air.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  pb.index <- pb.index + 1
  setTxtProgressBar(pb, pb.index)
  # ------------------------------------------
  # Modeling Temperature 
  # ------------------------------------------
  {
    # Load the saved model
    # Load the meta info for the betas
    
    betas.air_temperature <- ncdf4::nc_open(file.path(path.model, "air_temperature", "betas_air_temperature_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.air_temperature, "1"))
    ncdf4::nc_close(betas.air_temperature)
    
    dat.sim[["air_temperature"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      dat.temp$air_temperature = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==min(dat.mod$time.day)){
        sim.lag <- stack(lags.init$air_temperature)
        names(sim.lag) <- c("lag.air_temperature", "ens")
        
        sim.lag$lag.air_temperature_min <- stack(lags.init$air_temperature_min)[,1]
        sim.lag$lag.air_temperature_max <- stack(lags.init$air_temperature_max)[,1]
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["air_temperature"]][dat.mod$time.day==(i-1)  & dat.mod$hour==23,], dim=c(1, ncol(dat.sim$air_temperature)))))
        names(sim.lag) <- c("lag.air_temperature", "ens")
        sim.lag$lag.air_temperature_min <- stack(apply(dat.sim[["air_temperature"]][dat.mod$time.day==(i-1),], 2, min))[,1]
        sim.lag$lag.air_temperature_max <- stack(apply(dat.sim[["air_temperature"]][dat.mod$time.day==(i-1),], 2, max))[,1]
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # dat.temp$surface_downwelling_shortwave_flux_in_air <- stack(dat.sim$surface_downwelling_shortwave_flux_in_air[rows.now,])[,1]
      
      # Loading the saved model file
      load(file.path(path.model, "air_temperature", paste0("model_air_temperature_", day.now, ".Rdata")))
      mod.air_temperature.doy <- mod.save
      rm(mod.save)
      
      # Pull coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.air_temperature <- ncdf4::nc_open(file.path(path.model, "air_temperature", paste0("betas_air_temperature_", day.now, ".nc")))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.air_temperature, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.air_temperature)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.air_temperature.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$air_temperature), replace=T)
      
      for(j in 1:ncol(dat.sim$air_temperature)){
        dat.prop <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
        air_temperature_max.ens <- max(dat.temp[dat.temp$ens==paste0("X", j), "air_temperature_max.day"])
        air_temperature_min.ens <- min(dat.temp[dat.temp$ens==paste0("X", j), "air_temperature_min.day"])
        
        # Hard-coding in some bounds so we don't drift too far away from our given maxes & mins
        # Not going to worry for the moment about what happens if we undershoot out max/min since
        # it looks like we're normally decently close
        dat.prop[dat.prop>air_temperature_max.ens+2] <- air_temperature_max.ens+2
        dat.prop[dat.prop<air_temperature_min.ens-2] <- air_temperature_min.ens-2
        
        dat.sim[["air_temperature"]][rows.now,j] <- dat.prop
      }
      rm(mod.air_temperature.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  pb.index <- pb.index + 1
  setTxtProgressBar(pb, pb.index)
  
  # ------------------------------------------
  # Modeling precipitation.flux 
  # NOTE: For precipitation.flux, we're doing this differently:
  #       We're basicallly modeling the proportion of that day's precipitation as a function of 
  #       the hour of day.  Because our beta fitting method leverages the covariance among betas,
  #       we'll end up with a daily sum that's pretty close to our daily total (90-110%)
  # ------------------------------------------
  {
    # Load the meta info for the betas
    betas.precipitation.flux <- ncdf4::nc_open(file.path(path.model, "precipitation.flux", "betas_precipitation.flux_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.precipitation.flux, "1"))
    ncdf4::nc_close(betas.precipitation.flux)
    
    dat.sim[["precipitation.flux"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      dat.temp$precipitation.flux = 99999 # Dummy value so there's a column
      dat.temp$rain.prop = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==min(dat.mod$time.day)){
        sim.lag <- stack(lags.init$precipitation.flux)
        names(sim.lag) <- c("lag.precipitation.flux", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["precipitation.flux"]][dat.mod$time.day==(i-1)  & dat.mod$hour==23,], dim=c(1, ncol(dat.sim$precipitation.flux)))))
        names(sim.lag) <- c("lag.precipitation.flux", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      
      # Load the saved model
      load(file.path(path.model, "precipitation.flux", paste0("model_precipitation.flux_", day.now, ".Rdata")))
      mod.precipitation.flux.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.precipitation.flux <- ncdf4::nc_open(file.path(path.model, "precipitation.flux", paste0("betas_precipitation.flux_", day.now, ".nc")))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.precipitation.flux, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.precipitation.flux)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.precipitation.flux.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Re-distribute negative probabilities -- add randomly to make more peaky
      if(max(dat.pred)>0){ # If there's no rain on this day, skip the re-proportioning
        tmp <- 1:nrow(dat.pred) # A dummy vector of the 
        for(j in 1:ncol(dat.pred)){
          if(min(dat.pred[,j])>=0) next
          rows.neg <- which(dat.pred[,j]<0)
          rows.add <- sample(tmp[!tmp %in% rows.neg],length(rows.neg), replace=T)
          
          for(z in 1:length(rows.neg)){
            dat.pred[rows.add[z],j] <- dat.pred[rows.add[z],j] - dat.pred[rows.neg[z],j]
            dat.pred[rows.neg[z],j] <- 0  
          }
        }
        dat.pred <- dat.pred/rowSums(dat.pred)
        dat.pred[is.na(dat.pred)] <- 0
      }
      # Convert precip into real units
      dat.pred <- dat.pred*as.vector((dat.temp$precipitation.flux.day))
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$precipitation.flux), replace=T)
      
      for(j in 1:ncol(dat.sim$precipitation.flux)){
        dat.sim[["precipitation.flux"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      rm(mod.precipitation.flux.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  pb.index <- pb.index + 1
  setTxtProgressBar(pb, pb.index)
  
  # ------------------------------------------
  # Modeling surface_downwelling_longwave_flux_in_air 
  # ------------------------------------------
  {
    # Load the meta info for the betas
    betas.surface_downwelling_longwave_flux_in_air <- ncdf4::nc_open(file.path(path.model,"surface_downwelling_longwave_flux_in_air", "betas_surface_downwelling_longwave_flux_in_air_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.surface_downwelling_longwave_flux_in_air, "1"))
    ncdf4::nc_close(n.beta)
    
    dat.sim[["surface_downwelling_longwave_flux_in_air"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      dat.temp$surface_downwelling_longwave_flux_in_air = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==min(dat.mod$time.day)){
        sim.lag <- stack(lags.init$surface_downwelling_longwave_flux_in_air)
        names(sim.lag) <- c("lag.surface_downwelling_longwave_flux_in_air", "ens")
        
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["surface_downwelling_longwave_flux_in_air"]][dat.mod$time.day==(i-1)  & dat.mod$hour==23,], dim=c(1, ncol(dat.sim$surface_downwelling_longwave_flux_in_air)))))
        names(sim.lag) <- c("lag.surface_downwelling_longwave_flux_in_air", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # dat.temp$surface_downwelling_shortwave_flux_in_air <- stack(dat.sim$surface_downwelling_shortwave_flux_in_air[rows.now,])[,1]
      
      
      # Load the saved model
      load(file.path(path.model, "surface_downwelling_longwave_flux_in_air", paste0("model_surface_downwelling_longwave_flux_in_air_", day.now, ".Rdata")))
      mod.surface_downwelling_longwave_flux_in_air.doy <- mod.save
      rm(mod.list)
      
      
      # Pull the coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.surface_downwelling_longwave_flux_in_air <- ncdf4::nc_open(file.path(path.model,"surface_downwelling_longwave_flux_in_air", "betas_surface_downwelling_longwave_flux_in_air_1.nc"))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.surface_downwelling_longwave_flux_in_air, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.surface_downwelling_longwave_flux_in_air)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.surface_downwelling_longwave_flux_in_air.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      dat.pred <- dat.pred^2 # because squared to prevent negative numbers
      
      # Hard-coding some sanity bounds by ball-parking things from NLDAS & CRUNCEP
      # This is necessary if you have poorly constrained training models
      dat.pred[dat.pred<100] <- 100
      dat.pred[dat.pred>600] <- 600
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$surface_downwelling_longwave_flux_in_air), replace=T)
      
      for(j in 1:ncol(dat.sim$surface_downwelling_longwave_flux_in_air)){
        # test <- which(dat.temp$ens==paste0("X", j))
        dat.sim[["surface_downwelling_longwave_flux_in_air"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
      rm(mod.surface_downwelling_longwave_flux_in_air.doy) # Clear out the model to save memory
    }
    
  }
  # ------------------------------------------
  pb.index <- pb.index + 1
  setTxtProgressBar(pb, pb.index)
  # ------------------------------------------
  # Modeling air_pressure 
  # ------------------------------------------
  {
    # Load the meta info for the betas
    betas.air_pressure <- ncdf4::nc_open(file.path(path.model, "air_pressure", "betas_air_pressure_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.air_pressure, "1"))
    ncdf4::nc_close(betas.air_pressure)
    
    dat.sim[["air_pressure"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      dat.temp$air_pressure = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==min(dat.mod$time.day)){
        sim.lag <- stack(lags.init$air_pressure)
        names(sim.lag) <- c("lag.air_pressure", "ens")
        
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["air_pressure"]][dat.mod$time.day==(i-1)  & dat.mod$hour==23,], dim=c(1, ncol(dat.sim$air_pressure)))))
        names(sim.lag) <- c("lag.air_pressure", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # Load the saved model
      load(file.path(path.model, "air_pressure", paste0("model_air_pressure", day.now, ".Rdata")))
      mod.air_pressure.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.air_pressure <- ncdf4::nc_open(file.path(path.model, "air_pressure", paste0("betas_air_pressure_", day.now, ".nc")))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.air_pressure, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.air_pressure)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.air_pressure.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$air_pressure), replace=T)
      
      for(j in 1:ncol(dat.sim$air_pressure)){
        dat.sim[["air_pressure"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      rm(mod.air_pressure.doy) # Clear out the model to save memory
    }
  }

  # ------------------------------------------
  pb.index <- pb.index + 1
  setTxtProgressBar(pb, pb.index)
  # ------------------------------------------
  # Modeling specific_humidity 
  # ------------------------------------------
  {
    # Load the meta info for the betas
    betas.specific_humidity <- ncdf4::nc_open(file.path(path.model, "specific_humidity", "betas_specific_humidity_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.specific_humidity, "1"))
    ncdf4::nc_close(betas.specific_humidity)
    
    dat.sim[["specific_humidity"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      dat.temp$specific_humidity = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==min(dat.mod$time.day)){
        sim.lag <- stack(lags.init$specific_humidity)
        names(sim.lag) <- c("lag.specific_humidity", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["specific_humidity"]][dat.mod$time.day==(i-1)  & dat.mod$hour==23,], dim=c(1, ncol(dat.sim$specific_humidity)))))
        names(sim.lag) <- c("lag.specific_humidity", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # Load the saved model
      load(file.path(path.model, "specific_humidity", paste0("model_specific_humidity_", day.now, ".Rdata")))
      mod.specific_humidity.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.specific_humidity <- ncdf4::nc_open(file.path(path.model, "specific_humidity", paste0("betas_specific_humidity_", day.now, ".nc")))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.specific_humidity, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.specific_humidity)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.specific_humidity.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      dat.pred <- exp(dat.pred) # because log-transformed
      
      # specific_humidity sometimes ends up with high or infinite values, so lets make sure those get brought down a bit
      if(max(dat.pred)>0.03)  {
        specific_humidity.fix <- ifelse(quantile(dat.pred, 0.99)<0.03, quantile(dat.pred, 0.99), 0.03)
        dat.pred[dat.pred>specific_humidity.fix] <- specific_humidity.fix
      }
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$specific_humidity), replace=T)
      
      for(j in 1:ncol(dat.sim$specific_humidity)){
        dat.sim[["specific_humidity"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      rm(mod.specific_humidity.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  pb.index <- pb.index + 1
  setTxtProgressBar(pb, pb.index)
  # ------------------------------------------
  # Modeling wind_speed 
  # ------------------------------------------
  {
    # Load the meta info for the betas
    betas.wind_speed <- ncdf4::nc_open(file.path(path.model, "wind_speed", "betas_wind_speed_1.nc"))
    n.beta <- nrow(ncdf4::ncvar_get(betas.wind_speed, "1"))
    ncdf4::nc_close(betas.wind_speed)
    
    dat.sim[["wind_speed"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in min(dat.mod$time.day):max(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "air_temperature_max.day", "air_temperature_min.day", "precipitation.flux.day", "surface_downwelling_shortwave_flux_in_air.day", "surface_downwelling_longwave_flux_in_air.day", "air_pressure.day", "specific_humidity.day", "wind_speed.day",
                                     "next.air_temperature_max", "next.air_temperature_min", "next.precipitation.flux", "next.surface_downwelling_shortwave_flux_in_air", "next.surface_downwelling_longwave_flux_in_air", "next.air_pressure", "next.specific_humidity", "next.wind_speed")]
      dat.temp$wind_speed = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==min(dat.mod$time.day)){
        sim.lag <- stack(lags.init$wind_speed)
        names(sim.lag) <- c("lag.wind_speed", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["wind_speed"]][dat.mod$time.day==(i-1)  & dat.mod$hour==23,], dim=c(1, ncol(dat.sim$wind_speed)))))
        names(sim.lag) <- c("lag.wind_speed", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # Load the saved model
      load(file.path(path.model, "wind_speed", paste0("model_wind_speed_", day.now, ".Rdata")))
      mod.wind_speed.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      rows.beta <- sample(1:n.beta, n.ens, replace=T)
      betas.wind_speed <- ncdf4::nc_open(file.path(path.model, "wind_speed", paste0("betas_wind_speed_", day.now, ".nc")))
      Rbeta <- as.matrix(ncdf4::ncvar_get(betas.wind_speed, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      ncdf4::nc_close(betas.wind_speed)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.wind_speed.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      #dat.pred <- dat.pred^2 # because square-rooted to prevent negative numbers
      
      # Hard code an upper-level sanity check on the wind_speed
      # 20 m/s = 45 mph; not hurricane strength, but plenty strong enough for most models
      # dat.pred[dat.pred>20] <- 20
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$wind_speed), replace=T)
      
      for(j in 1:ncol(dat.sim$wind_speed)){
        dat.sim[["wind_speed"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      rm(mod.wind_speed.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  return(dat.sim)
  
}