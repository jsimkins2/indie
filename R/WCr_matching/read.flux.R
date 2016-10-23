read.flux <- function(fluxdata,
                       overwrite=FALSE, verbose=FALSE, ...){  
  library(ncdf4)
  library(zoo)
  var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                               "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                               "specific_humidity","precipitation_flux"),
                   units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
  )
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
  