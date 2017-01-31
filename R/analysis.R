library(udunits2)
library(PEcAn.all)
years = 2013:2015
obs <- list()
for(i in seq_along(years)){
  obs[[i]] <- load.L2Ameriflux.cf(paste0("WCr_data/raw/US-WCr.",years[i],".nc"))$NEE
}
obs = unlist(obs)
obs[obs==-9999] = NA
NEEo <- ud.convert(obs,"ug m-2 s-1","g m-2 s-1")/12*10e6
obs_mean <- mean(obs, na.rm=TRUE)
sum((NEEo*12/10e6*1800), na.rm = TRUE)
#----------------------------------MACA-----------------------------------------------------------

load("~/Cheas/WillowCreek/WCr_output/WCr_flux/PEcAn_99000000009/ensemble.ts.analysis.99000000168.NEE.2013.2015.Rdata")
flux_NEE = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#flux_NEE = flux_NEE - abs(obs_mean)
plot(flux_NEE, type = "l", col = "blue")

load("WCr_output/MACA/PEcAn_9900000008/ensemble.ts.analysis.99000000199.NEE.2013.2015.Rdata")
maca_deb = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#maca_deb = maca_deb - abs(obs_mean)
lines(maca_deb, col = "red")

load("WCr_output/MACA_raw/rcp85_1/ensemble.ts.analysis.99000000190.NEE.2013.2015.Rdata")
maca_raw = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#maca_raw = maca_raw - abs(obs_mean)
maca_raw30 = vector()
for (i in seq_along(maca_raw)){
  x_rep = vector()
  x_rep = rep(maca_raw[i],48)
  maca_raw30 = append(maca_raw30,x_rep)
}
plot(maca_raw30, col = "yellow", type = "l")
lines(maca_deb, col = "blue")
lines(flux_NEE, col = "orange")

df.m = data.frame(flux_NEE,maca_raw30,maca_deb)
boxplot(df.m,col = c("orange","yellow","blue"), outline=FALSE, names = c("Observed Forcing", "Raw Model", "Debiased & Downscaled"))

plot(maca_deb[(12000-48*6):12000], col = "dodgerblue4", type = "l", ylim = range(-200:200), lwd = 3, ylab = "umolC m-2 s-1", xlab = "Half Hour Intervals")
lines(flux_NEE[(12000-48*6):12000], col = "aquamarine4", lwd = 3)
lines(maca_raw30[(12000-48*6):12000], col = "coral2", lwd = 3)
title(main = "September 1st-7th MACA-GFDL-ESM2M NEE")
legend("bottomleft", legend = c("Observed Forcing","Raw Model","Debiased & Downscaled"),
       col = c("aquamarine4","coral2","dodgerblue4"), lty = c(1,1,1), lwd = 3, cex = 2)

flux_daily = colMeans(matrix(flux_NEE, nrow=48))
maca_daily = colMeans(matrix(maca_deb, nrow=48))
df.daily = data.frame(flux_daily,maca_raw,maca_daily)
df.daily = df.daily
boxplot(df.daily,col = c("aquamarine4","coral2","dodgerblue4"), outline=FALSE, names = c("Observed Forcing", "Raw Model", "Debiased & Downscaled"),
        ylab = "umolC m-2 s-1")
title(main = "MACA-GFDL-ESM2M NEE")

flux_sum = sum(flux_NEE)

#--------------------------------------------GFDL-------------------------------------------------

load("~/Cheas/WillowCreek/WCr_output/WCr_flux/PEcAn_99000000009/ensemble.ts.analysis.99000000168.NEE.2013.2015.Rdata")
flux_NEE = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#flux_NEE = flux_NEE - abs(obs_mean)
plot(flux_NEE, type = "l", col = "blue")

load("WCr_output/GFDL_debias_dwnsc/ESM2M/ensemble.ts.analysis.99000000201.NEE.2013.2015.Rdata")
gfdl_deb = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#gfdl_deb = gfdl_deb - abs(obs_mean)
#gfdl_deb[gfdl_deb > .0000004] = 0
lines(gfdl_deb)

load("WCr_output/GFDL_raw/ESM2M/ensemble.ts.analysis.99000000205.NEE.2013.2015.Rdata")
gfdl_raw = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#gfdl_raw = gfdl_raw - abs(obs_mean)
#gfdl_raw[gfdl_raw > .0000004] = 0
gfdl_raw30 = vector()
for (i in seq_along(gfdl_raw)){
  x_rep = vector()
  x_rep = rep(gfdl_raw[i],6)
  gfdl_raw30 = append(gfdl_raw30,x_rep)
}

plot(gfdl_raw30, col = "yellow", type = "l")
lines(gfdl_deb, col = "blue")
lines(flux_NEE, col = "orange")

df.m = data.frame(flux_NEE,gfdl_raw30,gfdl_deb)
boxplot(df.m,col = c("aquamarine4","coral2","dodgerblue4"), outline=FALSE, names = c("Observed Forcing", "Raw Model", "Debiased & Downscaled"),ylab = "umolC m-2 s-1")
title(main = "GFDL NEE Flux Boxplot")
plot(gfdl_deb[(12000-48*6.5-5):(12000-24)], col = "dodgerblue4", type = "l", ylim = range(-200:150), lwd = 3, ylab = "umolC m-2 s-1", xlab = "Half Hour Intervals")
lines(flux_NEE[(12000-48*6):12000], col = "aquamarine4", lwd = 3)
lines(gfdl_raw30[(12000-48*6.5-5):(12000-24)], col = "coral2", lwd = 3)
title(main = "September 1st-7th GFDL-ESM2M NEE")
legend("topleft", legend = c("Observed Forcing","Raw Daily Model","Debiased & Downscaled"),
       col = c("aquamarine4","coral2","dodgerblue4"), lty = c(1,1,1), lwd = 3, cex = 2)

flux_6hourly = colMeans(matrix(flux_NEE, nrow=6))
gfdl_6hourly = colMeans(matrix(gfdl_deb, nrow=6))
df.6hourly = data.frame(flux_6hourly,gfdl_raw,gfdl_6hourly)
df.6hourly = df.6hourly
boxplot(df.6hourly,col = c("aquamarine4","coral2","dodgerblue4"), outline=FALSE, names = c("Observed Forcing", "Raw Model", "Debiased & Downscaled"),
        ylab = "umolC m-2 s-1")
title(main = "GFDL-ESM2M NEE")

#---------------------------------------Future--------------------------------------------------------
load("WCr_output/GFDL_future/PEcAn_99000000013/ensemble.ts.analysis.99000000172.NEE.2060.2064.Rdata")
gfdl_2060 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#gfdl_2060 = gfdl_2060 - abs(obs_mean)

load("WCr_output/GFDL_future/PEcAn_99000000014/ensemble.ts.analysis.99000000174.NEE.2090.2094.Rdata")
gfdl_2090 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#gfdl_2090 = gfdl_2090 - abs(obs_mean)

load("WCr_output/MACA_future/ESM2M_2060_100/ensemble.ts.analysis.99000000237.NEE.2060.2064.Rdata")
maca_2060 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#maca_2060 = maca_2060 - abs(obs_mean)

load("WCr_output/MACA_future/ESM2M_2090_100/ensemble.ts.analysis.99000000239.NEE.2090.2094.Rdata")
maca_2090 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
#maca_2090 = maca_2090 - abs(obs_mean)

df.future = data.frame(gfdl_2060,maca_2060,gfdl_2090,maca_2090)
df.future = df.future
boxplot(df.future, col = c("darkolivegreen1","darkorange","darkolivegreen1","darkorange"), 
        names = c("GFDL 2060-2064", "MACA 2060-2064", "GFDL 2090-2094", "MACA 2090-2094"),
        ylab = "umolC m-2 s-1",outline=FALSE)
title(main = "Future GFDL vs. Future MACA")

#---------------------------------------GFDL through the years-----------------------
gfdl_deb_g = gfdl_deb
gfdl_2060_g = gfdl_2060
gfdl_2090_g = gfdl_2090
boxplot(x = gfdl_deb_g,gfdl_2060_g,gfdl_2090_g, col = c("darkolivegreen1","darkolivegreen1","darkolivegreen1"), 
        names = c("2013-2015", "2060-2064", "2090-2094"),
        ylab = "umolC m-2 s-1",outline=FALSE)
title(main = "GFDL Time Progression")

maca_deb_g = maca_deb
maca_2060_g = maca_2060
maca_2090_g = maca_2090
boxplot(x = maca_deb_g,maca_2060_g,maca_2090_g, col = c("darkorange","darkorange","darkorange"), 
        names = c("2013-2015", "2060-2064", "2090-2094"),
        ylab = "umolC m-2 s-1",outline=FALSE)
title(main = "MACA Time Progression")


#--------------------------------------annual means--------------------------------------------------------
flux_NEE_sum = flux_NEE*12/10e6*1800
gfdl_deb_sum = gfdl_deb*12/10e6*1800
maca_deb_sum = maca_deb*12/10e6*1800
gfdl_raw30_sum = gfdl_raw30*12/10e6*1800
maca_raw30_sum = maca_raw30*12/10e6*1800
gfdl_2060_sum = gfdl_2060*12/10e6*1800
maca_2060_sum = maca_2060*12/10e6*1800
gfdl_2090_sum = gfdl_2090*12/10e6*1800
maca_2090_sum = maca_2090*12/10e6*1800

ann_flux = c(sum(flux_NEE_sum[1:17520]),sum(flux_NEE_sum[17521:35040]),sum(flux_NEE_sum[35041:52560]))
ann_gfdl_raw = c(sum(gfdl_raw30_sum[1:17520]),sum(gfdl_raw30_sum[17521:35040]),sum(gfdl_raw30_sum[35041:52560]))
ann_gfdl = c(sum(gfdl_deb_sum[1:17520]),sum(gfdl_deb_sum[17521:35040]),sum(gfdl_deb_sum[35041:52560]))
gf_ann = data.frame(ann_flux,ann_gfdl_raw,ann_gfdl)
boxplot(gf_ann, ylim = range(700:2100),col = c("aquamarine4","coral2","dodgerblue4"), outline=FALSE, names = c("Observed Forcing", "Raw Model", "Debiased & Downscaled"),
        ylab = "gC m-2 yr-1")
title(main = "GFDL-ESM2M NEE")

ann_flux = c(sum(flux_NEE_sum[1:17520]),sum(flux_NEE_sum[17521:35040]),sum(flux_NEE_sum[35041:52560]))
ann_maca_raw = c(sum(maca_raw30_sum[1:17520]),sum(maca_raw30_sum[17521:35040]),sum(maca_raw30_sum[35041:52560]))
ann_maca = c(sum(maca_deb_sum[1:17520]),sum(maca_deb_sum[17521:35040]),sum(maca_deb_sum[35041:52560]))
gf_maca = data.frame(ann_flux,ann_maca_raw,ann_maca)
boxplot(gf_maca, ylim = range(0:2100),col = c("aquamarine4","coral2","dodgerblue4"), outline=FALSE, names = c("Observed Forcing", "Raw Model", "Debiased & Downscaled"),
        ylab = "gC m-2 yr-1")
title(main = "MACA-ESM2M NEE")


ann_gfdl2060 = c(sum(gfdl_2060_sum[1:17520]),sum(gfdl_2060_sum[17521:35040]),sum(gfdl_2060_sum[35041:52560]),sum(gfdl_2060_sum[52561:70080]),sum(gfdl_2060_sum[70081:87600]))
ann_gfdl2090 = c(sum(gfdl_2090_sum[1:17520]),sum(gfdl_2090_sum[17521:35040]),sum(gfdl_2090_sum[35041:52560]),sum(gfdl_2090_sum[52561:70080]),sum(gfdl_2090_sum[70081:87600]))
ann_maca2060 = c(sum(maca_2060_sum[1:17520]),sum(maca_2060_sum[17521:35040]),sum(maca_2060_sum[35041:52560]),sum(maca_2060_sum[52561:70080]),sum(maca_2060_sum[70081:87600]))
ann_maca2090 = c(sum(maca_2090_sum[1:17520]),sum(maca_2090_sum[17521:35040]),sum(maca_2090_sum[35041:52560]),sum(maca_2090_sum[52561:70080]),sum(maca_2090_sum[70081:87600]))
future = data.frame(ann_gfdl2060,ann_maca2060,ann_gfdl2090,ann_maca2090)

boxplot(future, ylim = range(550:2400),col = c("darkolivegreen1","darkorange","darkolivegreen1","darkorange"), 
        names = c("GFDL 2060-2064", "MACA 2060-2064", "GFDL 2090-2094", "MACA 2090-2094"), outline = FALSE,
        ylab = "gC m-2 yr-1")
title(main = "GFDL vs. MACA")
