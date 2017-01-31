#problem 2.4
setwd("~/Wisconsin/2016/Environmental Biophysics/Homework 1")
met = read.csv('Problem_2_4.csv',header = TRUE)
met = data.frame(met)
omega = pi/12
t = c(1:24)
gamma = vector()
for (z in 1:length(t)){
  gamma[[z]] = 0.44 - 0.46*sin(omega*t[z] +0.9)+0.11*sin(2*omega*t[z]+0.9)
}

max_3rd = 97.5
min_4 = as.numeric(as.character(min(met$Tair..ºF.[1:24])))
max_4 = max(met$Tair..ºF.[1:24])
min_5 = 73.1
max_17th = 102.5
min_18 = min(met$Tair..ºF.[25:48])
min_19 = 63.8
max_18 = max(met$Tair..ºF.[25:48])

Temp_4 = vector()
for (z in 1:length(gamma)){
  if (z>0 & z<=5){
    Temp_4[[z]] = max_3rd*gamma[[z]]+min_4*(1-gamma[[z]])
  }
  if (z>5 & z<=14){
      Temp_4[[z]] = max_4*gamma[[z]]+min_4*(1-gamma[[z]])
      }
  if (z>14 & z<=24){
      Temp_4[[z]] = max_4*gamma[[z]]+min_5*(1-gamma[[z]])
      }
  }

Temp_18 = vector()
for (z in 1:length(gamma)){
  if (z>0 & z<=5){
    Temp_18[[z]] = max_17th*gamma[[z]]+min_18*(1-gamma[[z]])
  }
  if (z>5 & z<=14){
    Temp_18[[z]] = max_18*gamma[[z]]+min_18*(1-gamma[[z]])
  }
  if (z>14 & z<=24){
    Temp_18[[z]] = max_18*gamma[[z]]+min_19*(1-gamma[[z]])
  }
}

met_4 =  vector(met$Tair..ºF.[1:24])
time = c(1:24)
met_18 = met$Tair..ºF.[25:48]
a = c("Observed", "Modeled")
plot(time,met_4,type="l",col="red", ylab = "Temperature (C)", xlab = "Time of Day (Hr)")
lines(time,Temp_4,col="green")
legend(1, 95, legend=c("Observed", "Modeled"),
       col=c("red", "green"), lty=1:1, cex=0.8)
title(main = "Problem 2.4 July 4")

a = c("Observed", "Modeled")
plot(time,met_18,type="l",col="red", ylab = "Temperature (C)", xlab = "Time of Day (Hr)")
lines(time,Temp_18,col="green")
legend(1, 95, legend=c("Observed", "Modeled"),
       col=c("red", "green"), lty=1:1, cex=0.8)
title(main = "Problem 2.4 July 18")
