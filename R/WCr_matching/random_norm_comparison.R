#normal distribution for future time periods
# rnorm(1, mean = daily.flux$air_temperature[1], sd = sd(flux$air_temperature[1:48]))

tem.met = data.frame()

#Here is where we go from Daily Resolution to 6-Hourly Resolution by placing in NAs 
for (n in 1:length(var$CF.name)){
  for (x in seq(from=0, to=1460, by=4)){
    tem.met[x,n] = model[x/4,n]}}

colnames(tem.met) = var$CF.name

for (n in seq_along(var$CF.name)){
  for (x in seq(from=1, to = 1460, by = 4)){
    tem.met[x,n] = rnorm(1, mean = model$air_temperature[x/4], sd = sd(flux[(x*12-12+1):(x*12),n]))

  }}

spline.met = data.frame()
spline.met = na.spline(tem.met[,1:3])
spline.met = na.spline(tem.met[,6:10])
spline.met = data.frame(spline.met)

spline.met$surface_downwelling_longwave_flux_in_air = NA
spline.met$air_pressure = NA

spline.met[,4] = NA
spline.met[,5] = NA
spline.met[,6:10] = na.spline(tem.met[,6:10])

colnames(spline.met) = var$CF.name

ranorm = data.frame()

for (n in seq_along(var$CF.name)){
  for (x in seq_along(1:1460)){
    if (x < 30*4){lowday = 0
    highday = x+30*4}
    if (x > 30*4){lowday = x-30*4 
    highday = 365*4}
    if( x +30*4 < 365*4 & x-30*4 > 0){lowday = x-30*4 
    highday = x+30*4}
    ranorm[x,n] = rnorm(1, mean = spline.met[x,n], sd = sd(flux[(lowday*12):(highday*12),n]))
    
  }}

colnames(ranorm) = var$CF.name
