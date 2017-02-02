library(PEcAn.all)
library(ncdf4)
library(lubridate)
library(REddyProc)

#-------------------

download.Ameriflux("US-WCr",outfolder = "Christy_code/raw_ameriflux",start_date = "2006-01-01",end_date = "2008-12-31",overwrite = FALSE,verbose = FALSE)
met2CF.Ameriflux("Christy_code/raw_ameriflux/",in.prefix = "US-WCr",outfolder = "Christy_code/CF_ameriflux",start_date = "2006-01-01",end_date = "2008-12-31",overwrite = FALSE,verbose = FALSE)
metgapfill(in.path = "Christy_code/CF_ameriflux/",in.prefix = "US-WCr",outfolder = "Christy_code/gapfilled_ameriflux",start_date = "2008-01-01",end_date = "2008-12-31",overwrite = FALSE,verbose = FALSE)

download.MACA(outfolder = "Christy_code/maca",start_date = "2006-01-01",end_date = "2008-12-31",site_id = 1,lat.in = 45.80592667,lon.in = -90.07985917, overwrite = FALSE, verbose = FALSE)


substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

input_met = c("Christy_code/gapfilled_ameriflux/US-WCr.2006.nc", "Christy_code/gapfilled_ameriflux/US-WCr.2007.nc")# "Christy_code/gapfilled_ameriflux/US-WCr.2008.nc")

sub_str <- substrRight(input_met, 7)
year <- substr(sub_str, 1, 4)
year <- as.numeric(year)

var <- data.frame(CF.name = c("air_temperature", "air_temperature_max", "air_temperature_min", 
                              "surface_downwelling_longwave_flux_in_air", "air_pressure", "surface_downwelling_shortwave_flux_in_air", 
                              "eastward_wind", "northward_wind", "specific_humidity", "precipitation_flux"), 
                  units = c("Kelvin", "Kelvin", "Kelvin", "W/m2", "Pascal", "W/m2", "m/s", 
                            "m/s", "g/g", "kg/m2/s"))


dat_data = data.frame()
for (i in 1:length(input_met)){
  raw_train_data <- list()
  tem <- ncdf4::nc_open(input_met[i])
  dim <- tem$dim
  for (j in seq_along(var$CF.name)) {
    if (exists(as.character(var$CF.name[j]), tem$var) == FALSE) {
      raw_train_data[[j]] = NA
    } else {
      raw_train_data[[j]] = ncdf4::ncvar_get(tem, as.character(var$CF.name[j]))
    }
  }
  raw_train_data = data.frame(raw_train_data)
  colnames(raw_train_data) <- var$CF.name 
  assign(paste0("dat_data", i),raw_train_data)
  
}

lat_raw_train_data <- as.numeric(ncdf4::ncvar_get(tem, "latitude"))
lon_raw_train_data <- as.numeric(ncdf4::ncvar_get(tem, "longitude"))
ncdf4::nc_close(tem)

dat.train = rbind(dat_data[i],dat_data2)

for(i in 1:ncol(dat.train)){
  dat.train[is.na(dat.train[,i]), i] <- mean(dat.train[,i], na.rm = TRUE)
}
