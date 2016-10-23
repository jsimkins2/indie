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
