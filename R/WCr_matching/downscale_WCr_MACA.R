library(ncdf4)
library(zoo)
var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                             "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                             "specific_humidity","precipitation_flux"),
                 units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
)

 read.flux <- function(fluxdata,
                      overwrite=FALSE, verbose=FALSE, ...){  
  #substrRight <- function(x, n){
  # substr(x, nchar(x)-n+1, nchar(x))
  #}
  #sub_str <- substrRight(flux, 7)
  #year <- substr(sub_str,1, 4)
  fluxdata = 'US-WCr.2006.nc'
  verbose = FALSE
  
  #Read in the two datasets, with the dimensions of the fluxrce dataset being named flux.list
  flux <- list()
  flux.list <- list()
  tem <- nc_open(fluxdata)
  dim <- tem$dim
  for (j in seq_along(var$CF.name)){
    if (var$CF.name[j] == "air_temperature_max" || var$CF.name[j] == "air_temperature_min"){
      flux[[j]] = NA}
    else {
      flux[[j]] <- ncvar_get(tem,as.character(var$CF.name[j]))
      flux.list[[j]] <- ncvar_def(name=as.character(var$CF.name[j]), 
                                  units=as.character(var$units[j]), 
                                  dim=dim, 
                                  missval=-999, 
                                  verbose=verbose)
    }}
  lat_flux <- as.numeric(ncvar_get(tem,"latitude"))
  lon_flux <- as.numeric(ncvar_get(tem,"longitude"))
  #year <- as.numeric(year)
  nc_close(tem)
  
  #Create dataframes from the lists of data pulled from the source/model and give them column names 
  flux <- data.frame(flux)
  colnames(flux) <- var$CF.name
  flux$air_temperature_max = flux$air_temperature
  flux$air_temperature_min = flux$air_temperature
  
  
  return(flux)
  
 }
 flux = read.flux('US-WCr.2006.nc')
 
 
 read.model <- function(modeldata,
                        overwrite=FALSE, verbose=FALSE, ...){  
   
   model <- list()
   tow <- nc_open(modeldata)
   for (j in seq_along(var$CF.name)){
     model[[j]] <- ncvar_get(tow,as.character(var$CF.name[j]))
   }
   nc_close(tow)
   
   model <- data.frame(model)
   colnames(model) <- var$CF.name
   
   return(model)}
 model <- read.model('MACA.IPSL-CM5A-LR.rcp85.r1i1p1.2007.nc')
 
 upscale.flux <- function(flux,
                          overwrite=FALSE, verbose=FALSE, ...){  
   
   var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air",
                                "air_pressure","surface_downwelling_shortwave_flux_in_air","eastward_wind",
                                "northward_wind","specific_humidity","precipitation_flux"),
                    units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s"))
   
   row = nrow(flux)
   #there are 48 measurements per day in 
   step = row/365
   
   daily.flux = data.frame()
   for (n in seq_along(var$CF.name)){
     for (x in 1:365){
       daily.flux[x,n] <- mean(flux[(x*step-step+1):(x*step),n])}
   }
   
   colnames(daily.flux) = var$CF.name
   #adding air_temperature_max and air_temperature_min
   for (x in 1:365){
     daily.flux$air_temperature_max[x] <- max(flux$air_temperature[(x*step-step+1):(x*step)])
     daily.flux$air_temperature_min[x] <- min(flux$air_temperature[(x*step-step+1):(x*step)])}
   
   
   return(daily.flux)
 }
 daily.flux <- upscale.flux(flux)
 
 downscale.model <- function(model,
                             overwrite=FALSE, verbose=FALSE, ...){  
   
   var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                                "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                                "specific_humidity","precipitation_flux"),
                    units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
   )
   days = 1:365
   hh = 1:17520
   downscaled.model = data.frame()
   for (n in seq_along(var$CF.name)){
     for (x in seq_along(days)){
       if(is.na(daily.flux[x,n])== FALSE & is.na(model[x,n]) == FALSE){
         #setting window boundaries 
         if (x < 30){lowday = 0
         highday = x+30}
         if (x > 30){lowday = x-30 
         highday = 365}
         if( x +30 < 365 & x-30 > 0){lowday = x-30 
         highday = x+30}
         #finding the closest daily model value to the closest daily flux value +/- 30 days
         coln = which.min(abs(daily.flux[lowday:highday,n] - model[x,n]))
         m = lowday + coln
         #Grabbing the 30 min flux data from specified model day
         downscaled.model[(x*48-48+1):(x*48),n] = flux[(m*48-23):(m*48+24),n]}
       else {
         if(is.na(daily.flux[x,n])== TRUE | is.na(model[x,n]) == TRUE){
           downscaled.model[(x*48-48+1):(x*48),n] = NA}}
     }
   }
   
   colnames(downscaled.model) = var$CF.name
   return(downscaled.model)
 }
 downscaled.model <- downscale.model(model)
 

 
 #plotting, the red is the WCr Flux Data, the blue is the downscaled model data. 
 par(mfrow=c(3,2))
 plot(flux$air_temperature, type = "l", col = "red",xlab = "Time", ylab = "Temperature (K)")
 lines(downscaled.model$air_temperature, col = "blue")
 plot(flux$surface_downwelling_shortwave_flux_in_air, type = "l", col = "red",xlab = "Time", ylab = "Shortwave")
 lines(downscaled.model$surface_downwelling_shortwave_flux_in_air, col = "blue")
 plot(flux$eastward_wind, type = "l", col = "red",xlab = "Time", ylab = "Eastward Wind")
 lines(downscaled.model$eastward_wind, col = "blue")
 plot(flux$northward_wind, type = "l", col = "red",xlab = "Time", ylab = "Northward Wind")
 lines(downscaled.model$northward_wind, col = "blue")
 plot(flux$specific_humidity, type = "l", col = "red",xlab = "Time", ylab = "Specific Humidity")
 lines(downscaled.model$specific_humidity, col = "blue")
 plot(flux$precipitation_flux, type = "l", col = "red",xlab = "Time", ylab = "Precipitation Flux")
 lines(downscaled.model$precipitation_flux, col = "blue")
 mtext("Willow Creek 2006 Flux Data", side = 3, line = -2, outer = TRUE)
 