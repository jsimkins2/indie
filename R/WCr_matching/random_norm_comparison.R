#normal distribution for future time periods
# rnorm(1, mean = daily.flux$air_temperature[1], sd = sd(flux$air_temperature[1:48]))
library(ncdf4)
library(zoo)

var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                             "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                             "specific_humidity","precipitation_flux"),
                 units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
)

verbose = FALSE
#Reading in the file
source_met = 'US-WCr.2006.nc'
flux <- list()
flux.list <- list()
tem <- nc_open(source_met)
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

sou.met = data.frame(flux)
colnames(sou.met) = var$CF.name
tem.met = data.frame()

reso = 1460
step = 48/4
six.met = data.frame()
for (n in 1:length(var$CF.name)){
  for (x in 1:reso){
    six.met[x,n] <- mean(sou.met[(x*step-step+1):(x*step),n])}
}

colnames(six.met) = var$CF.name

reso = 365
step = 48
daily.met = data.frame()
for (n in seq_along(var$CF.name)){
  for (x in 1:reso){
    daily.met[x,n] <- mean(sou.met[(x*step-step+1):(x*step),n])}
}

colnames(daily.met) = var$CF.name
tem.met = data.frame()

#Here is where we go from Daily Resolution to 6-Hourly Resolution by placing in NAs 
for (n in seq_along(var$CF.name)){
  for (x in seq(from=0, to=1460, by=4)){
    tem.met[x,n] = daily.met[x/4,n]}}

colnames(tem.met) = var$CF.name

hr <- 6
dayhr <- 24
step <- dayhr/hr
yrval <- step*365
flux_res <- nrow(flux)
w <- flux_res/yrval
#Here is where we go from Daily Resolution to 6-Hourly Resolution by placing in NAs 
tem.met <- data.frame()
for (n in seq_along(var$CF.name)){
  for (x in seq(from=0, to=yrval, by=step)){
    tem.met[x,n] <- model[x/4,n]}}

colnames(tem.met) <- var$CF.name

s <- list()
for (i in seq_along(var$CF.name)){
  if(all(is.na(tem.met[[i]]))== FALSE){
    s[[i]] <-  na.spline(tem.met[[i]])}
  if(all(is.na(tem.met[[i]])) == TRUE){
    s[[i]] <- NA}}

spline.met <- data.frame(s)
colnames(spline.met) <- var$CF.name

rand.met <- data.frame()
for (n in seq_along(var$CF.name)){
  for (x in seq_along(1:yrval)){
    if(all(is.na(spline.met[n]))== TRUE){
      rand.met[x,n] <- NA}
    if(all(is.na(spline.met[n]))== FALSE){
      if (x < wd*step){lowday <- 0
      highday <- x+wd*step}
      if (x > wd*step){lowday <- x-wd*step
      highday <- 365*step}
      if( x +wd*step < 365*step & x-wd*4 > 0){lowday <- x-wd*step
      highday <- x+wd*step}
      rand.met[x,n] <- rnorm(1, mean = spline.met[x,n], sd = .5*sd(spline.met[(lowday*w):(highday*w),n]))}
    
  }}

colnames(rand.met) <- var$CF.name
#Increase the deviation of Surface Downwelling
for (x in seq_along(1:1460)){
  if (x < wd*4){lowday <- 0
  highday <- x+wd*4}
  if (x > wd*4){lowday <- x-wd*4 
  highday <- 365*4}
  if( x +wd*4 < 365*4 & x-wd*4 > 0){lowday <- x-wd*4 
  highday <- x+wd*4}
  rand.met$surface_downwelling_shortwave_flux_in_air[x] <- rnorm(1, mean = spline.met$surface_downwelling_shortwave_flux_in_air[x],
                                                                 sd = .5*sd(flux$surface_downwelling_shortwave_flux_in_air[(lowday*w):(highday*w)]))}

rand.met$surface_downwelling_shortwave_flux_in_air[rand.met$surface_downwelling_shortwave_flux_in_air < 0] <- 0
rand.met$precipitation_flux[rand.met$precipitation_flux < 0] <- 0
rand.met$specific_humidity[rand.met$specific_humidity < 0] <- 0
