
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}



heat_deb = gfdl_deb[22:52560]
heat_deb = append(heat_deb,rep(0,23))
gf = data.frame(c(1:48))
step = 48
for (i in seq_len(365*3)) {
  n <- heat_deb[(i*step - step+1):(step*i)]
  gf[,i] <- n
}

gf[gf>200]=226.3889
gf = as.matrix(gf)


gf = t(gf)
pal.1=colorRampPalette(c("blue", "white", "red"), space="rgb")

png("gf1.png", width=6, height=6, units="in", res=200)
layout(matrix(c(1,2,3,0,4,0), nrow=3, ncol=3), widths=c(4,4,1), heights=c(4,1))
layout.show(4)

#1st image
breaks <- seq(min(gf), max(gf),length.out=100)
par(mar=c(1,1,1,1))
image(seq(dim(gf)[1]), seq(dim(gf)[2]), gf, 
      col=pal.1(length(breaks)-1), breaks=breaks, xaxt="n" , yaxt="n", ylab="Time of Day", xlab="Day of Year")

par(mar=c(3,1,1,1))
image.scale(gf, col=pal.1(length(breaks)-1), breaks=breaks, horiz=TRUE)
box()


dev.off()


library(lattice) 
levelplot(gf[1:ncol(gf),ncol(gf):20])
