
# This is just the names of the variables in the file
maca = data.frame(DAP.name = c("air_temperature", "rsds","uas","vas","huss","pr"),
                  CF.name = c("air_temperature","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"),
                  units = c('Kelvin',"W/m2","m/s","m/s","g/g","kg/m2/s")
)

verbose = FALSE
#Reading in the file
source_met = 'US-WCr.2006.nc'
sou = list()
sou.list = list()
tem = nc_open(source_met)
dim = tem$dim
for (j in 1:length(maca$CF.name)){
  sou[[j]] <- ncvar_get(tem,as.character(maca$CF.name[j]))
  sou.list[[j]] <- ncvar_def(name=as.character(maca$CF.name[j]), units=as.character(maca$units[j]), dim=dim, missval=-999, verbose=verbose)
}

sou.met = data.frame(sou)
colnames(sou.met) = maca$CF.name
tem.met = data.frame()

reso = 1460
step = 48/4
six.met = data.frame()
for (n in 1:length(maca$DAP.name)){
  for (x in 1:reso){
    six.met[x,n] <- mean(sou.met[(x*step-step+1):(x*step),n])}
}
colnames(six.met) = c("air_temperature","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux")

reso = 365
step = 48
daily.met = data.frame()
for (n in 1:length(maca$DAP.name)){
  for (x in 1:reso){
    daily.met[x,n] <- mean(sou.met[(x*step-step+1):(x*step),n])}
}

colnames(daily.met) = c("air_temperature","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux")

#Here is where we go from Daily Resolution to 6-Hourly Resolution by placing in NAs 
tem.met = data.frame()
for (n in 1:length(maca$DAP.name)){
  for (x in seq(from=0, to=1460, by=4)){
    tem.met[x,n] = daily.met[x/4,n]}}

colnames(tem.met) = maca$CF.name

# Spline fitting the NAs for every variabile but Precipitation
spline.met = data.frame()
spline.met = na.spline(tem.met[,1:5])
spline.met = data.frame(spline.met)
#Now we randomly assign a time for the Precipitation data from the original file to be placed
r = list()
r = sample(0:4,365,replace=T)
for (x in seq(from=4, to=1460, by=4)){
  spline.met[x-r[x/4],6] = tem.met[x/4,6]}
colnames(spline.met) = maca$CF.name
spline.met$precipitation_flux[is.na(spline.met$precipitation_flux)] <- 0

'''
par(mfrow=c(3,2))
plot(c(1:1460),six.met$air_temperature, type = "l", col = "red",xlab = "Time", ylab = "Temperature (K)")
lines(spline.met$air_temperature, col = "blue")
plot(c(1:1460),six.met$surface_downwelling_shortwave_flux_in_air, type = "l", col = "red",xlab = "Time", ylab = "Shortwave")
lines(spline.met$surface_downwelling_shortwave_flux_in_air, col = "blue")
plot(c(1:1460),six.met$eastward_wind, type = "l", col = "red",xlab = "Time", ylab = "Eastward Wind")
lines(spline.met$eastward_wind, col = "blue")
plot(c(1:1460),six.met$northward_wind, type = "l", col = "red",xlab = "Time", ylab = "Northward Wind")
lines(spline.met$northward_wind, col = "blue")
plot(c(1:1460),six.met$specific_humidity, type = "l", col = "red",xlab = "Time", ylab = "Specific Humidity")
lines(spline.met$specific_humidity, col = "blue")
plot(c(1:1460),six.met$precipitation_flux, type = "l", col = "red",xlab = "Time", ylab = "Precipitation Flux")
lines(spline.met$precipitation_flux, col = "blue")
mtext("Willow Creek 2006 Flux Data", side = 3, line = -2, outer = TRUE)
'''

par(mfrow=c(1,1))
plot(c(1:28),six.met$air_temperature[1:28], type = "l", col = "red", xlab = "Time", ylab = "Temperature (K)")
lines(spline.met$air_temperature[1:28], col = "blue")
legend("topright", legend =  c("Spline Interpolation","6-hour Means", "Daily Flux Value"), cex =0.8, 
       col=c("blue","red", "orange"), lty=c(1,1,NA), pch = c(NA, NA, 1))
points(x = tem.met$air_temperature[1:28], col = "orange")
title(main = "1 Week of Willow Creek Data")


 #p = paste(as.character(v[1:730]), collapse=",")
     # v = (trunc(six.met$surface_downwelling_shortwave_flux_in_air, prec = 1))
      