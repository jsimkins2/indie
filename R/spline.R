
# This is just the names of the variables in the file
maca = data.frame(DAP.name = c("air_temperature", "rsds","uas","vas","huss","pr"),
                  CF.name = c("air_temperature","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"),
                  units = c('Kelvin',"W/m2","m/s","m/s","g/g","kg/m2/s")
)

verbose = FALSE
#Reading in the file
source_met = 'MACA.IPSL-CM5A-LR.rcp85.r1i1p1.2007.nc'
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

#Here is where we go from Daily Resolution to 6-Hourly Resolution by placing in NAs 
for (n in 1:length(maca$DAP.name)){
  for (x in seq(from=0, to=1460, by=4)){
    tem.met[x,n] = sou.met[x/4,n]}}

colnames(tem.met) = maca$CF.name

# Spline fitting the NAs for every variabile but Precipitation
spline.met = data.frame()
spline.met = na.spline(tem.met[,1:5])
spline.met = data.frame(spline.met)
#Now we randomly assign a time for the Precipitation data from the original file to be placed
r = list()
r = sample(0:4,365,replace=T)
for (x in seq(from=4, to=1460, by=4)){
  spline.met[x-r[x/4],6] = sou.met[x/4,6]}
colnames(spline.met) = maca$CF.name
spline.met$precipitation_flux[is.na(spline.met$precipitation_flux)] <- 0


