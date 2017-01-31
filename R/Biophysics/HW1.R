#James Simkins
#Homework 1
#2.1A

setwd("~/Wisconsin/EnvBio")


height = vector()
height = c(6.4, 3.2, 1.6, 0.8, 0.4, 0.2, 0.1)
mean_air = vector()
mean_air = c(31.08, 31.72, 32.37, 33.05, 33.80, 34.74, 36.91)


plot(mean_air, height, main = "2.1 A",xlab = "Mean Air Temp (C)", ylab = "Height (m)")

#2.1B
ln.vector = vector()
for (z in 1:length(height)){
  ln.vector[[z]]= log((height[[z]]-(0.6*.15))/(0.02*.15))
}

plot(mean_air, ln.list, main = "2.1 B", xlab = "Mean Air Temp (C)", ylab = "ln[(z-d)/zh)")

data = data.frame(height, mean_air, ln.vector)

#2.3
depth = c(6, 12, 18, 24, 30)

tz = vector()
for (d in 1:length(depth)){
  tz[[d]]=25 - 20^(-depth[[d]]/10)
}
tmin = vector()
for (d in 1:length(depth)){
  tmin[[d]]=25 - 20^(-depth[[d]]/10)
}

plot(depth, tz, main = "2.3", xlab = "depth (m)", ylab = "Min Temperature (C)")



#2.4

met = read.csv('Problem_2_4.csv', header = TRUE)
