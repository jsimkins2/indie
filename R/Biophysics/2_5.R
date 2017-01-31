# Problem 2.5
setwd("~/Wisconsin/2016/Environmental Biophysics/Homework 1")
data = read.csv('2_5_data.csv')
df = data.frame(data)

GDD.2015 = vector()
tb= 9.1
tbase = vector()
for (j in 1:length(df$Tmax.2015)){
  if (is.na(df$Tmax.2015[j]) == FALSE && is.na(df$Tmin.2015[j]) == FALSE){
      if (df$Tmax.2015[j] > 30){
        GDD.2015[[j]] = 0}
      if (df$Tmax.2015[j] < 30){
        if (((df$Tmax.2015[j] + df$Tmin.2015[j])/2) < 9.1){
          tbase[[j]] = (df$Tmax.2015[j] + df$Tmin.2015[j])/2}
        if (((df$Tmax.2015[j] + df$Tmin.2015[j])/2) > 9.1){
          tbase[[j]] = 9.1
        }
        GDD.2015[[j]] = (df$Tmax.2015[j] + df$Tmin.2015[j])/2 - tbase[j]
      }
    
  }
  if (is.na(df$Tmax.2015[j]) == TRUE || is.na(df$Tmin.2015[j]) == TRUE){
    GDD.2015[[j]] <- NA
  }
}

s.2015 = vector()
for (j in 1:length(df$Tmax.2015)){
  if (is.na(GDD.2015[j]) == TRUE){
    GDD.2015[j] = 0}
  s.2015[[j]] = sum(GDD.2015[1:j])
}
#method 2 
day = c(1:365)

plot(day,s.2010[1:365],type="l",col="red", ylab = "Growing Degree Days", xlab = "Day Of Year", ylim = range(0:4000))
lines(day,s.2011[1:365],col="green")
lines(day,s.2012[1:365],col="black")
lines(day,s.2013[1:365],col="pink")
lines(day,s.2014[1:365],col="blue")
lines(day,s.2015[1:365],col="orange")

legend("topleft", legend=c("2015", "2015", "2015", "2015", "2015","2015","Leaf Emergence","Silking","Physiological Maturity"),
       col=c("red", "green","black","pink","blue","orange","gray60","yellow","brown"), lty=1:1, cex=0.8)
title(main = "Method 2")
abline(h = 61, col = "gray60")
abline(h = 650,col = "yellow")
abline(h = 1300, col = "brown")
