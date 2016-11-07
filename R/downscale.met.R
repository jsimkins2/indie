##' downscale.met takes model data and flux data from the same site and temporally 
##'    downscales the modeled data based on statistics from the flux data.
##' @name downscale.met
##' @title downscale.met
##' @export
##' @param outfolder
##' @param modeldata - the modeled dataset in NC format. Could be CMIP5 dataset (see download.GFDL or download.MACA) or other modeled dataset 
##' @param fluxdata - the observed dataset in NC format i.e. Flux Tower dataset (see download.Fluxnet2015 or download.Ameriflux) 
##' @param method - how you want to downscale. Default will be 'normal', which performs a random normal distribution on the 
##'                 spline interpolated points using the value as the mean and the standard deviation of the flux data +/- window_days from the value 
##' @param site.id
##' @param ensemble_members - a numeric value that selects how many ensemble members you would like to run 
##' @param window_days - a numeric value setting the size of the window for standard deviation calculations
##' @param resolution - resolution in HOURS that you would like the modeled data to be downscaled to, default option assumes 6 hourly, i.e. a vector is created 
##'                     with values every 6 hours  
##' @param swdn_method - Downscaling downwelling shortwave flux in air presents an interesting problem. Through our attempts of solving it, we found 2 methods. 
##'                      The 'sine' method is our standard, and fits a spline interpolation over a daily average value during the daylight length. The 'Waichler'
##'                      method comes from Waichler and Wigtosa 2003 and presents a method on calculating surface downwelling shortwave flux in air based on
##'                      temperature mins/maxes and precipitation. 
##' @author James Simkins
downscale.met <- function(outfolder, modeldata, fluxdata, site_id, method='normal',ensemble_members = 10, window_days = 20,
                          resolution = 6, swdn_method = "sine", overwrite=FALSE, verbose=FALSE, ...){  
  reso <- resolution
  wd <- window_days
  w <- 24/reso
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  sub_str <- substrRight(modeldata, 7)
  year <- substr(sub_str,1, 4)
  
  # eph_year for shortwave calculations. Since MACA doesn't contain leap years, eph_year is set to 2006 is year given is a leap year. 
  eph_year <- year
  mod_name <- substr(modeldata,1,4)
  
  
  fluxdata <- "US-WCr.2006.nc"
  #Variable names 
  var <- data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                               "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                               "specific_humidity","precipitation_flux"),
                   units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
  )
  #Reading in the flux data
    flux <- list()
    tem <- ncdf4::nc_open(fluxdata)
    dim <- tem$dim
    for (j in seq_along(var$CF.name)){
      if (exists(as.character(var$CF.name[j]),tem$var) == FALSE){
        flux[[j]] = NA}
      else {
        flux[[j]] <- ncdf4::ncvar_get(tem,as.character(var$CF.name[j]))}
    }
    lat_flux <- as.numeric(ncdf4::ncvar_get(tem,"latitude"))
    lon_flux <- as.numeric(ncdf4::ncvar_get(tem,"longitude"))
    ncdf4::nc_close(tem)
    
    flux <- data.frame(flux)
    colnames(flux) <- var$CF.name
    flux$air_temperature_max <- flux$air_temperature
    flux$air_temperature_min <- flux$air_temperature
  
  #Reading in the model data
    model <- list()
    tow <- ncdf4::nc_open(modeldata)
    for (j in seq_along(var$CF.name)){
      if (exists(as.character(var$CF.name[j]),tow$var) == FALSE){
        model[[j]] = NA}
      else {
        model[[j]] <- ncdf4::ncvar_get(tow,as.character(var$CF.name[j]))}
      }
    ncdf4::nc_close(tow)
    
    model <- data.frame(model)
    colnames(model) <- var$CF.name
  
  #Going to upscale the observations to the desired downscale resolution so we can take the standard deviatiion of the residuals
    # MACA doesn't have leap days, so this parses the flux data by removing the leap day
    if (substr(modeldata,1,4) == 'MACA'){
    if (length(flux$air_temperature)%%365 > 0){
      floor(length(flux$air_temperature)/365)}
      eph_year = 2006}
    
    reso_len <- 365*24/reso
    step <- length(flux$air_temperature)/reso_len
    upscale_flux <- data.frame()
    for (n in 1:length(var$CF.name)){
      for (x in 1:reso_len){
        upscale_flux[x,n] <- mean(flux[(x*step-step+1):(x*step),n])}
    }
    colnames(upscale_flux) <- var$CF.name
    for (x in 1:reso_len){
      upscale_flux$air_temperature_max[x] <- max(flux$air_temperature[(x*step-step+1):(x*step)])
      upscale_flux$air_temperature_min[x] <- min(flux$air_temperature[(x*step-step+1):(x*step)])}
    
    # Now we start a for loop for the ensemble members and begin downscaling. If method = "normal", a random normal distribution
    # is used to downscale as so; (mean = value of modeled data) (sd = +/- window_days of upscale_flux data at the same time intervals)
    # Future work is to implement more methods for downscaling
    
    for (e in seq_along(1:ensemble_members)){
      if (method == "normal"){
        # Ephemeris is a function to calculate sunrise/sunset times and daylength for SW calculations
        ephemeris <- function(lat, lon, date, span=1, tz="UTC") {
          
          lon.lat <- matrix(c(lon, lat), nrow=1)
          
          # using noon gets us around daylight saving time issues
          day <- as.POSIXct(sprintf("%s 12:00:00", date), tz=tz)
          sequence <- seq(from=day, length.out=span , by="days")
          
          sunrise <- maptools::sunriset(lon.lat, sequence, direction="sunrise", POSIXct.out=TRUE)
          sunset <- maptools::sunriset(lon.lat, sequence, direction="sunset", POSIXct.out=TRUE)
          solar_noon <- maptools::solarnoon(lon.lat, sequence, POSIXct.out=TRUE)
          
          data.frame(date=as.Date(sunrise$time),
                     sunrise=as.numeric(format(sunrise$time, "%H%M")),
                     solarnoon=as.numeric(format(solar_noon$time, "%H%M")),
                     sunset=as.numeric(format(sunset$time, "%H%M")),
                     day_length=as.numeric(sunset$time-sunrise$time)) 
          
        }
        
      #### Temperature 
      randtemp <- vector()
      for (x in seq_along(model$air_temperature)){
        four <- vector()
        for (n in seq_along(1:4)){
          lowday <- x*w+n-wd
          highday <- lowday+wd*2+n
          if (lowday < 0){lowday = n}
          if (highday > reso_len){highday = reso_len-4+n}
          four[n] <- rnorm(1, mean = model$air_temperature[x],
                            sd = sd(upscale_flux$air_temperature[seq(from=lowday,to=highday,by=4)]))}
        randtemp <- append(randtemp,four)}
      
      ### Air Temperature Max and Min
      temp_max <- vector()
      for (x in seq_along(model$air_temperature)){
        four <- vector()
        for (n in seq_along(1:4)){
          lowday <- x*w+n-wd
          highday <- lowday+wd*2+n
          if (lowday < 0){lowday = n}
          if (highday > reso_len){highday = reso_len-4+n}
          four[n] = rnorm(1, mean = model$air_temperature_max[x],
                          sd = sd(upscale_flux$air_temperature_max[seq(from=lowday,to=highday,by=4)]))}
        temp_max <- append(temp_max,four)}
      
      temp_min <- vector()
      for (x in seq_along(model$air_temperature)){
        four <- vector()
        for (n in seq_along(1:4)){
          lowday <- x*w+n-wd
          highday <- lowday+wd*2+n
          if (lowday < 0){lowday = n}
          if (highday > reso_len){highday = reso_len-4+n}
          four[n] = rnorm(1, mean = model$air_temperature_min[x],
                          sd = sd(upscale_flux$air_temperature_min[seq(from=lowday,to=highday,by=4)]))}
        temp_min <- append(temp_min,four)}
      #### Precipitation
      # this takes the daily total of precipitation and uses that as a total possible amount of precip. 
      # It randomly distributes the values of precipitation 
      rand_vect_cont <- function(N, M, sd = 1) {
        vec <- rnorm(N, M/N, sd)
        vec / sum(vec) * M
      }
      precip <- vector()
      for (x in seq_along(model$precipitation_flux)){
        lowday <- (x-wd)*w
        highday <- (x+wd)*w
        if (lowday < 0){lowday = 0}
        if (highday > reso_len){highday = reso_len}
        four <- vector()
        four <- rand_vect_cont(4,model$precipitation_flux[x], 
                               sd = sd(upscale_flux$precipitation_flux[lowday:highday]))
        four[four<0] = 0
        precip <- append(precip,four)}
      
      #### Specific Humidity
      spechum <- vector()
      for (x in seq_along(model$specific_humidity)){
        lowday <- (x-wd)*w
        highday <- (x+wd)*w
        if (lowday < 0){lowday = 0}
        if (highday > reso_len){highday = reso_len}
        four <- vector()
        for (n in seq_along(1:4)){
          four[n] <- rnorm(1, mean = model$specific_humidity[x], 
                          sd = sd(upscale_flux$specific_humidity[lowday:highday]))}
        spechum <- append(spechum,four)}
      spechum[spechum < 0] = 0
      
      #### Winds
      
      east <- vector()
      for (x in seq_along(model$eastward_wind)){
        lowday <- (x-wd)*w
        highday <- (x+wd)*w
        if (lowday < 0){lowday = 0}
        if (highday > reso_len){highday = reso_len}
        four <- vector()
        for (n in seq_along(1:4)){
          four[n] <- rnorm(1, mean = model$eastward_wind[x], 
                          sd = sd(upscale_flux$eastward_wind[lowday:highday]))}
        east <- append(east,four)}
  
      north <- vector()
      for (x in seq_along(model$northward_wind)){
        lowday <- (x-wd)*w
        highday <- (x+wd)*w
        if (lowday < 0){lowday = 0}
        if (highday > reso_len){highday = reso_len}
        four <- vector()
        for (n in seq_along(1:4)){
          four[n] <- rnorm(1, mean = model$northward_wind[x],
                          sd = sd(upscale_flux$northward_wind[lowday:highday]))}
        north <- append(north,four)}
  
      #### Shortwave Downwelling Radiation
      swmodel <- model$surface_downwelling_shortwave_flux_in_air
      swdn <- vector()
      # The sine method produces an hourly sine wave of 
      if (lubridate::leap_year(eph_year) == TRUE){
        span = 366}
      if (lubridate::leap_year(eph_year)){
        span = 365}
      if (swdn_method == "sine"){
        eph <- ephemeris(lat_flux,lon_flux,date = paste0(eph_year,"-01-01 00:00:00"), 
                       span = span, tz = "UTC")
        day_len <- eph$day_length
      for (i in seq_along(swmodel)){
        t <- seq(from=pi/day_len[i],to=pi,by=pi/day_len[i])
        wav <- ((swmodel[i]*(24/day_len[i]))/0.637)*sin(t)
        
        srs <- eph$sunrise
        
        hr <- substr(srs[i],1,2)
        hr <-as.numeric(hr)
        #6 is for CST time zone which is UTC-6. Need to automatically find TZ by latitude.
        hr <- hr-6
        
        l <- vector()
        for (n in seq_along(1:hr)){
          l[n] = 0}
        for (n in seq_along(wav)){
          l[n+hr] = wav[n]}
        for(n in seq_along(1:(24-(length(wav)+hr)))){
          l[n+hr+length(wav)] = 0}
        
        swdn <- append(swdn,l)}
        
      swflux <- vector()
      sw_step <- length(swdn)/reso_len
      for (x in 1:reso_len){
        swflux[x] <- mean(swdn[(x*sw_step-sw_step+1):(x*sw_step)])}}
      # The Waichler method doesn't need averaged SW flux values, it models SW downwelling flux based on Tmax-Tmin and Precipitation
      # Reference is Waichler and Wigtosa 2003. Our no-precip coefficient is 2 instead of 1 (fits the data better)
      if (swdn_method == "Waichler"){
        inter <- paste0(reso," hour")
        days <- seq(as.POSIXct("2006-01-01 00:00:00", tz = "US/Central"), as.POSIXct("2006-12-31 18:00:00", tz = "US/Central"), by=inter)
        
        Z <- Ratmosphere::SZA(days,lat_flux, lon_flux)
        I <- 1000*cos_d(Z)
        m <- vector()
        for (i in seq_along(1:12)){
          m[i] <- Hmisc::monthDays(as.Date(paste0(year,'-',i,'-01')))}
        bmlist <- vector()
  
        Bm <- c(0.2089,0.2857, 0.2689, 0.2137, 0.1925, 0.2209, 0.2527, 0.2495,0.2232, 0.1728, 0.1424, 0.1422)
        for (x in seq_along(Bm)){
          mlen <- list()
          mlen <- rep(Bm[x],m[x]*4)
          bmlist <- append(bmlist,mlen)}
        A <-.73
        C <- .7
        hdry <- vector()
        for (i in seq_along(downscale.met$air_temperature)){
          if (downscale.met$precipitation_flux[i] > 0){p=.65}
          if (downscale.met$precipitation_flux[i] == 0){p=2}
          hdry[i] <- A*p*(1-exp(-1*bmlist[i]*((model$air_temperature_max[i]-model$air_temperature_min[i])^C)))
        }
        hdry[hdry<0]=0
        swflux <- hdry*I
        swflux[swflux<0]=0
        }
      #### Longwave Downwelling 
      if (all(is.na(model$surface_downwelling_longwave_flux_in_air)) == TRUE){
        lwflux<-vector()
        lwflux<- rep(NA,reso_len)}
      if (all(is.na(model$surface_downwelling_longwave_flux_in_air)) == FALSE){
      lwflux<- vector()
      for (x in seq_along(model$surface_downwelling_longwave_flux_in_air)){
        lowday <- (x-wd)*w
        highday <- (x+wd)*w
        if (lowday < 0){lowday = 0}
        if (highday > reso_len){highday = reso_len}
        four <- vector()
        for (n in seq_along(1:4)){
          four[n] <- rnorm(1, mean = model$surface_downwelling_longwave_flux_in_air[x], 
                          sd = sd(upscale_flux$surface_downwelling_longwave_flux_in_air[lowday:highday]))}
        lwflux <- append(lwflux,four)}}
      
      #### Atmospheric Pressure
      
      if (all(is.na(model$air_pressure)) == TRUE){
        pres<-vector()
        pres<- rep(NA,reso_len)}
      if (all(is.na(model$air_pressure)) == FALSE){
        pres<- vector()
        for (x in seq_along(model$air_pressure)){
          lowday <- (x-wd)*w
          highday <- (x+wd)*w
          if (lowday < 0){lowday = 0}
          if (highday > reso_len){highday = reso_len}
          four <- vector()
          for (n in seq_along(1:4)){
            four[n] <- rnorm(1, mean = model$air_pressure[x], sd = sd(upscale_flux$air_pressure[lowday:highday]))}
          pres <- append(pres,four)}}
      
      } #this ends the method = 'normal' downscaling phase

    #### Putting all the variables together in a data frame
    downscaled.met <- data.frame(randtemp,temp_max,temp_min,lwflux,pres,swflux,east,north,spechum,precip)
    colnames(downscaled.met) <- var$CF.name

    flux.list <- list()
    lat <- ncdf4::ncdim_def(name='latitude', units='degree_north', vals=lat_flux, create_dimvar=TRUE)
    lon <- ncdf4::ncdim_def(name='longitude', units='degree_east', vals=lon_flux, create_dimvar=TRUE)
    time <- ncdf4::ncdim_def(name='time', units="sec", vals=(1:reso_len)*reso*3600, create_dimvar=TRUE, unlim=TRUE)
    dim<-list(lat,lon,time)
    
    for (j in seq_along(var$CF.name)){
    flux.list[[j]] <- ncdf4::ncvar_def(name=as.character(var$CF.name[j]), 
                                units=as.character(var$units[j]), 
                                dim=dim, 
                                missval=-999, 
                                verbose=verbose)}
    
    rows <- 1
  dir.create(outfolder, showWarnings=FALSE, recursive=TRUE)
  results <- data.frame(file=character(rows), host=character(rows),
                        mimetype=character(rows), formatname=character(rows),
                        startdate=character(rows), enddate=character(rows),
                        dbfile.name = paste0(mod_name,".dwnsc.ens",e,".",year,sep="."),
                        stringsAsFactors = FALSE)
  

  loc.file <- file.path(outfolder, paste0(mod_name,".dwnsc.ens",e,".",year,".nc"))
  
  loc <- ncdf4::nc_create(filename=loc.file, vars=flux.list, verbose=verbose)
  for(j in seq_along(var$CF.name)){
    ncdf4::ncvar_put(nc=loc, varid=as.character(var$CF.name[j]), vals=downscaled.met[[j]])
  }
  ncdf4::nc_close(loc)
  
  results$file <- loc.file
  results$host <- PEcAn.utils::fqdn()
  results$startdate <- paste0(year,"-01-01 00:00:00", tz = "UTC")
  results$enddate <- paste0(year,"-12-31 23:59:59", tz = "UTC")
  results$mimetype <- 'application/x-netcdf'
  results$formatname <- 'CF Meteorology'
  
  invisible(results)
    }
}

# downscale.met('trial', 'MACA.IPSL-CM5A-LR.rcp85.r1i1p1.2006.nc', 'US-WCr.2006.nc')
  