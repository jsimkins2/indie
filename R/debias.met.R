##' debias.met takes source_met, sample_met and output new_met that debiases source_met against sample_met
##' @name debias.met
##' @title debias.met
##' @export
##' @param outfolder
##' @param source_met - the met data that you want to be altered 
##' @param sample_met - the met data that we use as the trainer(the more accurate dataset)
##' @param site.id
##' @author James Simkins
debias.met <- function(outfolder, source_met, sample_met, site_id, overwrite=FALSE, verbose=FALSE, ...){  
  require(PEcAn.utils)
  require(lubridate)
  require(ncdf4)
  outfolder = paste0(outfolder,"_site_",paste0(site_id %/% 1000000000, "-", site_id %% 1000000000))
  
  var = data.frame(DAP.name = c("tas","rlds","ps","rsds","uas","vas","huss","pr"),
                   CF.name = c("air_temperature","surface_downwelling_longwave_flux_in_air","air_pressure","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"),
                   units = c('Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
  )
  
  #Need to pull the year from the string filename
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  sub_str= substrRight(source_met, 7)
  year = substr(sub_str,1, 4)
  #Read in the two datasets, with the dimensions of the source dataset being named sou.list
  sou = list()
  sou.list = list()
  tem = nc_open(source_met)
  dim = tem$dim
  for (j in 1:length(var$CF.name)){
    sou[[j]] = ncvar_get(tem,as.character(var$CF.name[j]))
    sou.list[[j]] = ncvar_def(name=as.character(var$CF.name[j]), units=as.character(var$units[j]), dim=dim, missval=-999, verbose=verbose)
  }
  nc_close(tem)
  
  sam = list()
  tow = nc_open(sample_met)
  for (j in 1:length(var$CF.name)){
    sam[[j]] = ncvar_get(tow,as.character(var$CF.name[j]))
  }
  nc_close(tow)
  
  #Create dataframes from the lists of data pulled from the source/sample and give them column names 
  sou = data.frame(sou)
  sam = data.frame(sam)
  colnames(sou) = c("tas","rlds","ps","rsds","uas","vas","huss","pr")
  colnames(sam) = c("tas","rlds","ps","rsds","uas","vas","huss","pr")
  
  #Grab the means of the source and sample, find the difference, and correct the source dataset accordingly
  mean_sou = apply(sou,2,mean)
  mean_sam = apply(sam,2,mean)
  mean_diff = mean_sam - mean_sou
  debi = list()
  for (k in 1:length(mean_diff)){
    debi[[k]] = (sou[[k]] + mean_diff[[k]])
  }
  
  rows = 1
  dir.create(outfolder, showWarnings=FALSE, recursive=TRUE)
  results <- data.frame(file=character(rows), host=character(rows),
                        mimetype=character(rows), formatname=character(rows),
                        startdate=character(rows), enddate=character(rows),
                        dbfile.name = paste("debias.met",sep="."),#"GFDL",
                        stringsAsFactors = FALSE)
  
  debi = data.frame(debi)
  loc.file = file.path(outfolder,paste("debias",year,"nc",sep="."))
  
  loc <- nc_create(filename=loc.file, vars=sou.list, verbose=verbose)
  for(j in 1:nrow(var)){
    ncvar_put(nc=loc, varid=as.character(var$CF.name[j]), vals=debi[[j]])
  }
  nc_close(loc)
  
  results$file <- loc.file
  results$host <- fqdn()
  results$startdate <- paste0(year,"-01-01 00:00:00")
  results$enddate <- paste0(year,"-12-31 23:59:59")
  results$mimetype <- 'application/x-netcdf'
  results$formatname <- 'CF Meteorology'
  
  invisible(results)
}


