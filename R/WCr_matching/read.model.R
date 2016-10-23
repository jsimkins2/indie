read.model <- function(modeldata,
           overwrite=FALSE, verbose=FALSE, ...){  
library(ncdf4)
library(zoo)
var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                             "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                             "specific_humidity","precipitation_flux"),
                 units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
)
modeldata = 'MACA.IPSL-CM5A-LR.rcp85.r1i1p1.2007.nc'

model <- list()
tow <- nc_open(modeldata)
for (j in seq_along(var$CF.name)){
  model[[j]] <- ncvar_get(tow,as.character(var$CF.name[j]))
}
nc_close(tow)

model <- data.frame(model)
colnames(model) <- var$CF.name

return(model)}