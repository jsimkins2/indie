require(ggplot2)
require(ncdf4)
require(lubridate)
require(data.table)
setwd("~/Wisconsin/Thesis/Debias")

var = data.frame(DAP.name = c("tas","rlds","ps","rsds","uas","vas","huss","pr"),
                 CF.name = c("air_temperature","surface_downwelling_longwave_flux_in_air","air_pressure","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"),
                 units = c('Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
)

gfdl = list()
tem = nc_open('GFDL.CM3.rcp45.r1i1p1.2006.nc')
for (j in 1:length(var$CF.name)){
  gfdl[[j]] = ncvar_get(tem,as.character(var$CF.name[j]))
}
tower = list()
tow = nc_open('US-WCr.2006.nc')
for (j in 1:length(var$CF.name)){
  tower[[j]] = ncvar_get(tow,as.character(var$CF.name[j]))
}

gfdl = data.frame(gfdl)
tower = data.frame(tower)
colnames(gfdl) = c("tas","rlds","ps","rsds","uas","vas","huss","pr")
colnames(tower) = c("tas","rlds","ps","rsds","uas","vas","huss","pr")

gdates = seq(from = as.POSIXct("2006-01-01 00:00"), by = "3 hours", length.out = 2920)
gdates.names <- strftime(gdates, "%Y-%m-%d %H:%M:%S", tz="GMT")
row.names(gfdl) <- gdates.names
fdates = seq(from = as.POSIXct("2006-01-01 00:00:00"), by = "30 min", length.out = nrow(tower))
fdates.names <- strftime(fdates, "%Y-%m-%d %H:%M:%S", tz="GMT")
row.names(tower) <- fdates.names

mean_t = apply(tower,2,mean)
mean_g = apply(gfdl,2,mean)

mean_diff = mean_t - mean_g

test = list()
for (k in 1:length(mean_diff)){
  test[[k]] = (gfdl[[k]] + mean_diff[[k]])
}

gfdl = data.frame(test)
colnames(gfdl) = c("tas","rlds","ps","rsds","uas","vas","huss","pr")
gdates = seq(from = as.POSIXct("2006-01-01 00:00"), by = "3 hours", length.out = 2920)
gdates.names <- strftime(gdates, "%Y-%m-%d %H:%M:%S", tz="GMT")
row.names(gfdl) <- gdates.names

