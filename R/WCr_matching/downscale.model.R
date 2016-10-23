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
