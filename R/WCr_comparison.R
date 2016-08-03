require(ggplot2)
require(ncdf4)

setwd("~/Wisconsin/Thesis/Data")
var = data.frame(DAP.name = c("tas","rlds","ps","rsds","uas","vas","huss","pr"),
                 CF.name = c("air_temperature","surface_downwelling_longwave_flux_in_air","air_pressure","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"),
                 units = c('Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
)

gfdl = list()
tem = nc_open('GFDL.CM3.rcp45.r1i1p1.2006.nc')
for (j in 1:length(var$CF.name)){
  gfdl[[j]] = ncvar_get(tem,as.character(var$CF.name[j]))
}
flux = read.csv('FLX_US-WCr_FLUXNET2015_SUBSET_HH_1999-2014_1-1.csv')

#Only looking at 2006
flux = flux[which(flux$TIMESTAMP_START > 200600000000 & flux$TIMESTAMP_START < 200700000000), ]

flux$TA_F = flux$TA_F + 273.15

gx = seq(1,2920,by=1)
plot(flux$TIMESTAMP_START,flux$TA_F,type="l",col="red", xlab = "Year",
     ylab = "Temperature (K)",main = "GFDL vs. Tower")
lines(gx,gfdl[[1]],col="green")
grid(NULL,NULL)
