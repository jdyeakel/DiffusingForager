#Starvation dynamics

library(deSolve)
#library(rgl)
library(RColorBrewer)
library(wesanderson)
source("R/filled_contour.r")
source("R/smooth_pal.R")

#Starving forager ODE
source("R/starv_forage_ode.R")

state <- c(
  R = 0.5,            #Resource
  H = 0.5,          #Starvers
  Full = 0.5)        #Non-Starvers

time <- seq(1,1000, by = 0.1)

parameters <- c(
  alpha = 0.5,
  lambda = 0.2,   
  sigma = 4,
  rho = 0.5,
  beta = 0.2,
  mu = 0.2,
  epsilon=0.001
)

out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
#Plot ALL
colors <- brewer.pal(5,"Set1")
plot(out[,1],out[,2], xlab = "time", ylab = "-",type="l",lwd=2,col=colors[2],ylim=c(0,max(as.numeric(cbind(out[,2],out[,3],out[,4]))))) #Blue
lines(out[,1],out[,3],col=colors[5],type="l",lwd=2) #Orange
lines(out[,1],out[,4],col=colors[3],type="l",lwd=2) #Green

#Extinction analysis across different sigma values
sigma_reps <- 20
sigma_seq <- seq(0.3,4,length.out=sigma_reps)
reps <- 10
extinct <- matrix(0,reps,length(sigma_seq))
for (r in 1:reps) {
  for (i in 1:sigma_reps) {
    state <- c(
      R = 0.5,            #Resource
      H = 0.5,          #Starvers
      Full = 0.5)        #Non-Starvers
    
    time <- seq(1,1000, by = 0.1)
    
    parameters <- c(
      alpha = 0.5,
      lambda = 0.2,   
      sigma = sigma_seq[i],
      rho = 0.5,
      beta = 0.2,
      mu = 0.2,
      epsilon=0.001
    )
    out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
    threshold <- 0.2
    if (any(out[,3]+out[,4] < threshold)) {
      extinct[r,i] <- 1
    } else {
      extinct[r,i] <- 0
    }
    print(paste("Finished r=",r," i=",i))
  }
}
ext_pr <- apply(extinct,2,sum)/reps
plot(sigma_seq,ext_pr,type="b",pch=16,xlab="Starvation rate",ylab="Extinction probability")
