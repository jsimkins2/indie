load("~/Cheas/WillowCreek/WCr_output/WCr_flux/PEcAn_99000000009/ensemble.ts.analysis.99000000168.NEE.2013.2015.Rdata")
flux_NEE = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

load("WCr_output/MACA/PEcAn_9900000008/ensemble.ts.analysis.99000000199.NEE.2013.2015.Rdata")
maca_deb = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

load("WCr_output/MACA_raw/rcp85_1/ensemble.ts.analysis.99000000190.NEE.2013.2015.Rdata")
maca_raw = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
maca_raw30 = vector()
for (i in seq_along(maca_raw)){
  x_rep = vector()
  x_rep = rep(maca_raw[i],48)
  maca_raw30 = append(maca_raw30,x_rep)
}

load("WCr_output/GFDL_debias_dwnsc/ESM2M/ensemble.ts.analysis.99000000201.NEE.2013.2015.Rdata")
gfdl_deb = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

load("WCr_output/GFDL_raw/ESM2M/ensemble.ts.analysis.99000000205.NEE.2013.2015.Rdata")
gfdl_raw = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6
gfdl_raw30 = vector()
for (i in seq_along(gfdl_raw)){
  x_rep = vector()
  x_rep = rep(gfdl_raw[i],6)
  gfdl_raw30 = append(gfdl_raw30,x_rep)
}

load("WCr_output/GFDL_future/ESM2M_2060_100/ensemble.ts.analysis.99000000247.NEE.2060.2064.Rdata")
gfdl_2060 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

load("WCr_output/GFDL_future/ESM2M_2090_100/ensemble.ts.analysis.99000000249.NEE.2090.2094.Rdata")
gfdl_2090 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

load("WCr_output/MACA_future/ESM2M_2060_100/ensemble.ts.analysis.99000000237.NEE.2060.2064.Rdata")
maca_2060 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

load("WCr_output/MACA_future/ESM2M_2090_100/ensemble.ts.analysis.99000000239.NEE.2090.2094.Rdata")
maca_2090 = ud.convert(ensemble.ts.analysis$mean, "kg m-2 s-1","g m-2 s-1")/12*10e6

df.future = data.frame(gfdl_2060,maca_2060,gfdl_2090,maca_2090)

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

ann_maca_raw = c(sum(maca_raw30_sum[1:17520]),sum(maca_raw30_sum[17521:35040]),sum(maca_raw30_sum[35041:52560]))
ann_maca = c(sum(maca_deb_sum[1:17520]),sum(maca_deb_sum[17521:35040]),sum(maca_deb_sum[35041:52560]))
maca_ann = data.frame(ann_flux,ann_maca_raw,ann_maca)

ann_gfdl2060 = c(sum(gfdl_2060_sum[1:17520]),sum(gfdl_2060_sum[17521:35040]),sum(gfdl_2060_sum[35041:52560]),sum(gfdl_2060_sum[52561:70080]),sum(gfdl_2060_sum[70081:87600]))
ann_gfdl2090 = c(sum(gfdl_2090_sum[1:17520]),sum(gfdl_2090_sum[17521:35040]),sum(gfdl_2090_sum[35041:52560]),sum(gfdl_2090_sum[52561:70080]),sum(gfdl_2090_sum[70081:87600]))
ann_maca2060 = c(sum(maca_2060_sum[1:17520]),sum(maca_2060_sum[17521:35040]),sum(maca_2060_sum[35041:52560]),sum(maca_2060_sum[52561:70080]),sum(maca_2060_sum[70081:87600]))
ann_maca2090 = c(sum(maca_2090_sum[1:17520]),sum(maca_2090_sum[17521:35040]),sum(maca_2090_sum[35041:52560]),sum(maca_2090_sum[52561:70080]),sum(maca_2090_sum[70081:87600]))
future_sum = data.frame(ann_gfdl2060,ann_maca2060,ann_gfdl2090,ann_maca2090)

#---------------------departure plotting

m1 = mean((ann_gfdl2060 - 1846.341) *.1)
m2 = mean((ann_maca2060 - 1846.341) *.1)
m3 = mean((ann_gfdl2090 - 1846.341) *.1)
m4 = mean((ann_maca2090 - 1846.341) *.1)

m5 = mean((ann_gfdl_raw - 1846.341)*.1)
m6 = mean((ann_gfdl - 1846.341)*.1)
m7 = mean((ann_maca_raw - 1846.341)*.1)
m8 = mean((ann_maca - 1846.341)*.1)

fut = list(m1,m2,m3,m4)
fut = unlist(fut)
png("future_gfdl_maca.png", width=10, height=6,res = 200, units = "in")
#par(mar=c(2,2,2,2))
bp = barplot(fut, ylim = range(-150:100),col = c("dodgerblue4","coral1","dodgerblue4","coral1"), 
             names = c("GFDL 2060-2064", "MACA 2060-2064", "GFDL 2090-2094", "MACA 2090-2094"),
             ylab = "gC m-2 yr-1", cex.axis = 1.4, cex.names = 1.4, cex.lab = 1.4)
text(x=bp[1], y=fut[1], labels=paste0("+",round(fut[1],0)," gC m-2 yr-1"), pos=3, xpd=NA, cex = 1.4)
text(x=bp[2], y=fut[2], labels=paste0(round(fut[2],0)," gC m-2 yr-1"), pos=1, xpd=NA, cex = 1.4)
text(x=bp[3], y=fut[3], labels=paste0("+",round(fut[3],0)," gC m-2 yr-1"), pos=3, xpd=NA, cex = 1.4)
text(x=bp[4], y=fut[4], labels=paste0(round(fut[4],0)," gC m-2 yr-1"), pos=1, xpd=NA, cex = 1.4)
text(x=bp[1], y=fut[1], labels="30% increase", pos=1, xpd=NA, cex = 1.4)
mtext(paste0((mean(ann_gfdl)-mean(ann_gfdl2060))/mean(ann_gfdl),"%"),side = 1,line = 2,at = .75)
abline(a = 0,b = 0)
dev.off()

pres = list(m5,m6,m7,m8)
pres = unlist(pres)

png("pres_gfdl_maca.png", width=10, height=6,res = 200, units = "in")
bp = barplot(pres, ylim = range(-200:25),col = c("dodgerblue4","dodgerblue4","coral1","coral1"), 
             names = c("GFDL Native", "Debiased/Downscaled", "MACA Native", "Debiased/Downscaled"),
             ylab = "gC m-2 yr-1", cex.axis = 1.4, cex.names = 1.4, cex.lab = 1.4)
abline(a = 0,b = 0)
text(x=bp, y=pres, labels=paste0(round(pres,0)," gC m-2 yr-1"), pos=1, xpd=NA, cex = 1.4)
dev.off()

#-----------------now doing future bsed off of changes from gfdl regular

m1 = (mean(ann_gfdl2060) - mean(ann_gfdl)) *.1
m2 = (mean(ann_maca2060) - mean(ann_maca)) *.1
m3 = (mean(ann_gfdl2090) - mean(ann_gfdl)) *.1
m4 = (mean(ann_maca2090) - mean(ann_maca)) *.1
fut = list(m1,m2,m3,m4)
fut = unlist(fut)
png("future_gfdl_maca.png", width=10, height=6,res = 200, units = "in")
#par(mar=c(2,2,2,2))
bp = barplot(fut, ylim = range(-100:200),col = c("dodgerblue4","coral1","dodgerblue4","coral1"), 
             names = c("GFDL 2060-2064", "MACA 2060-2064", "GFDL 2090-2094", "MACA 2090-2094"),
             ylab = "gC m-2 yr-1", cex.axis = 1.4, cex.names = 1.4, cex.lab = 1.4)
text(x=bp[1], y=fut[1], labels=paste0("+",round(fut[1],0)," gC m-2 yr-1"), pos=3, xpd=NA, cex = 1.4)
text(x=bp[2], y=fut[2], labels=paste0(round(fut[2],0)," gC m-2 yr-1"), pos=1, xpd=NA, cex = 1.4)
text(x=bp[3], y=fut[3], labels=paste0("+",round(fut[3],0)," gC m-2 yr-1"), pos=3, xpd=NA, cex = 1.4)
text(x=bp[4], y=fut[4], labels=paste0(round(fut[4],0)," gC m-2 yr-1"), pos=1, xpd=NA, cex = 1.4)
text(x=bp[1], y=fut[1], labels=paste0("+",as.integer(m1/(mean(ann_gfdl)*.1)*100),"%"), pos=1, xpd=NA, cex = 1.4)
text(x=bp[2], y=fut[2], labels=paste0(as.integer(m2/(mean(ann_maca)*.1)*100),"%"), pos=3, xpd=NA, cex = 1.4)
text(x=bp[3], y=fut[3], labels=paste0("+",as.integer(m3/(mean(ann_gfdl)*.1)*100),"%"), pos=1, xpd=NA, cex = 1.4)
text(x=bp[4], y=fut[4], labels=paste0(as.integer(m4/(mean(ann_maca)*.1)*100),"%"), pos=3, xpd=NA, cex = 1.4)

abline(a = 0,b = 0)
dev.off()


#---------------------------adding % to present plot 

png("pres_gfdl_maca.png", width=10, height=6,res = 200, units = "in")
bp = barplot(pres, ylim = range(-200:25),col = c("dodgerblue4","dodgerblue4","coral1","coral1"), 
             names = c("GFDL Native", "Debiased/Downscaled", "MACA Native", "Debiased/Downscaled"),
             ylab = "gC m-2 yr-1", cex.axis = 1.4, cex.names = 1.4, cex.lab = 1.4)
abline(a = 0,b = 0)
text(x=bp, y=pres, labels=paste0(round(pres,0)," gC m-2 yr-1"), pos=1, xpd=NA, cex = 1.4)
text(x = bp, y = pres, labels = paste0(as.integer(pres/(mean(ann_flux)*.1)*100),"%"), pos=3,xpd=NA,cex = 1.4)
dev.off()
