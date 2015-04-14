#Starvation dynamics

library(deSolve)
library(RColorBrewer)

#Simulate Grassland-Sapling-Tree-Woodland dynamics
source("R/starv_forage_ode.R")

pvec <- seq(0,5,0.1)
lpvec <- length(pvec)
ResSS <- numeric(lpvec)
ResSD <- numeric(lpvec)
PopSS <- numeric(lpvec)
PopSD <- numeric(lpvec)


for (i in 1:lpvec) {
  
state <- c(
  R = 0.5,            #Resource
  X = 0.5,          #Starvers
  Y = 0.5)        #Non-Starvers

parameters <- c(
  alpha = 1,     
  p = pvec[i],       
  mu_x = 0.5,     
  mu_y = 0.25,
  lambda = 8
)

time <- seq(0,500, by = 0.1)

out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)

ResSS[i] <- mean(out[4000:5000,2])
ResSD[i] <- sd(out[4000:5000,2])

PopSS[i] <- mean(out[4000:5000,3] + out[4000:5000,4])
PopSD[i] <- sd(out[4000:5000,3] + out[4000:5000,4])

}

plot(pvec,ResSS,ylim=c(0,1),type="b")
plot(pvec,ResSD,ylim=c(0,1),type="b")
plot(pvec,ResSD/ResSS,ylim=c(0,max(ResSD/ResSS)),type="b")

plot(pvec,PopSS,ylim=c(0,1),type="b")
plot(pvec,PopSD,ylim=c(0,1),type="b")
plot(pvec,PopSD/PopSS,ylim=c(0,max(PopSD/PopSS)),type="b")

#Plot ALL
colors <- brewer.pal(5,"Set1")
plot(out[,1],out[,2], xlab = "time", ylab = "-",type="l",lwd=2,col=colors[2],ylim=c(0,1)) #Blue
lines(out[,1],out[,3],col=colors[5],type="l",lwd=2) #Orange
lines(out[,1],out[,4],col=colors[3],type="l",lwd=2) #Green