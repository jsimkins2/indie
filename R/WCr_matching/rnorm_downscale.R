#normal distribution for future time periods
# rnorm(1, mean = daily.flux$air_temperature[1], sd = sd(flux$air_temperature[1:48]))
library(ncdf4)
library(zoo)

var = data.frame(CF.name = c("air_temperature","air_temperature_max","air_temperature_min","surface_downwelling_longwave_flux_in_air","air_pressure",
                             "surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind",
                             "specific_humidity","precipitation_flux"),
                 units = c('Kelvin','Kelvin','Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
)
yrval <- 1460
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

tem.met <- data.frame()
for (n in seq_along(var$CF.name)){
  for (x in seq(from=0, to=1460, by=4)){
    tem.met[x,n] <- daily.met[x/4,n]}}
colnames(tem.met) <- var$CF.name

###
####
#####
######### Temperature Flux
w <- 4
wd <- 10
randtemp = list()
for (x in seq_along(1:365)){
    lowday <- (x-wd)*w
    highday <- (x+wd)*w
    if (lowday < 0){lowday = 0}
    if (highday > 1460){highday = 1460}
  four <- list()
  for (n in seq_along(1:4)){
    four[[n]] = rnorm(1, mean = daily.met$air_temperature[x], sd = sd(six.met$air_temperature[lowday:highday]))
    st = sort(unlist(four))}
  randtemp = append(randtemp,st)}

randtemp = unlist(randtemp)

###
####
#####
########## shortwave 

w <- 4
wd <- 10
short = list()
for (x in seq_along(1:365)){
  lowday <- (x-wd)*w
  highday <- (x+wd)*w
  if (lowday < 0){lowday = 0}
  if (highday > 1460){highday = 1460}
  four <- list()
  for (n in seq_along(1:2)){
    four[[n]] = rnorm(1, mean = daily.met$surface_downwelling_shortwave_flux_in_air[x], sd = sd(six.met$surface_downwelling_shortwave_flux_in_air[lowday:highday]))
    four[[3]] = 0
    four[[4]] = 0
    four = unlist(four)
    four[four<0] = 0
    four[four > max(six.met$surface_downwelling_shortwave_flux_in_air[lowday:highday] + 10)] = max(six.met$surface_downwelling_shortwave_flux_in_air[lowday:highday] + 10)
    st = sort(four)}
  short = append(short,st)}

short = unlist(short)

###
####
#####
########## precipitation
# this takes the daily total of precipitation and uses that as a total possible amount of precip. 
# It randomly distributes the values of precipitation 
rand_vect_cont <- function(N, M, sd = 1) {
  vec <- rnorm(N, M/N, sd)
  vec / sum(vec) * M
}
w <- 4
wd <- 10
precip = list()
for (x in seq_along(1:365)){
  lowday <- (x-wd)*w
  highday <- (x+wd)*w
  if (lowday < 0){lowday = 0}
  if (highday > 1460){highday = 1460}
  four <- list()
  four <- rand_vect_cont(4,daily.met$precipitation_flux[x])
  four[four<0] = 0
  st = unlist(four)
  precip = append(precip,st)}

precip = unlist(precip)

###
####
#####
######### Specific Humidity

spechum = list()
for (x in seq_along(1:365)){
  lowday <- (x-wd)*w
  highday <- (x+wd)*w
  if (lowday < 0){lowday = 0}
  if (highday > 1460){highday = 1460}
  four <- list()
  for (n in seq_along(1:4)){
    four[[n]] = rnorm(1, mean = daily.met$specific_humidity[x], sd = sd(six.met$specific_humidity[lowday:highday]))
    st = unlist(four)
    st[st<0] = 0}
  spechum = append(spechum,st)}

spechum = unlist(spechum)


###
####
#####
######### Windy Winds
west = list()
for (x in seq_along(1:365)){
  lowday <- (x-wd)*w
  highday <- (x+wd)*w
  if (lowday < 0){lowday = 0}
  if (highday > 1460){highday = 1460}
  four <- list()
  for (n in seq_along(1:4)){
    four[[n]] = rnorm(1, mean = daily.met$eastward_wind[x], sd = sd(six.met$eastward_wind[lowday:highday]))
    st = unlist(four)}
  west = append(west,st)}

west = unlist(west)


north = list()
for (x in seq_along(1:365)){
  lowday <- (x-wd)*w
  highday <- (x+wd)*w
  if (lowday < 0){lowday = 0}
  if (highday > 1460){highday = 1460}
  four <- list()
  for (n in seq_along(1:4)){
    four[[n]] = rnorm(1, mean = daily.met$northward_wind[x], sd = sd(six.met$northward_wind[lowday:highday]))
    st = unlist(four)}
  north = append(north,st)}

north = unlist(north)

#Idea...we can run multiple random normal generators and increase the standard deviation with each decade for example to 
#get a better sense of how things will be in the future...this will 'increase the tails' of the distribution...
#we will also be able to use different distributions!!!
