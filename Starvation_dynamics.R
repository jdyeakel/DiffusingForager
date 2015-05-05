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
  X = 0.5,          #Starvers
  Y = 0.5)        #Non-Starvers

parameters <- c(
  alpha = 2, 
  epsilon = 0.5,
  sigma = 0.1,
  rho = 0.5,
  gamma = 2,    
  mu = 0.02
)

time <- seq(0,5000, by = 0.1)

out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
#Plot ALL
colors <- brewer.pal(5,"Set1")
plot(out[,1],out[,2], xlab = "time", ylab = "-",type="l",lwd=2,col=colors[2],ylim=c(0,max(as.numeric(cbind(out[,2],out[,3],out[,4]))))) #Blue
lines(out[,1],out[,3],col=colors[5],type="l",lwd=2) #Orange
lines(out[,1],out[,4],col=colors[3],type="l",lwd=2) #Green

#Plot R/X+Y vs. gammaC
with(as.list(parameters),{
  R <- out[,2]
  X <- out[,3]
  Y <- out[,4]
  C <- X + Y
  gammaC <- gamma*Y*R - 2*mu*X - mu*Y
  plot((R/C),gammaC,pch=".")
}
)




#Simulate dynamics over BOTH sigma and gamma
source("R/starv_forage_ode.R")

sigmavec <- seq(0,0.5,0.05)
l_sigmavec <- length(sigmavec)
gammavec <- seq(0.1,1,0.1)
l_gammavec <- length(gammavec)

ResSS_m <- matrix(0,l_sigmavec,l_gammavec)
ResSD_m <- matrix(0,l_sigmavec,l_gammavec)
XSS_m <- matrix(0,l_sigmavec,l_gammavec)
XSD_m <- matrix(0,l_sigmavec,l_gammavec)
YSS_m <- matrix(0,l_sigmavec,l_gammavec)
YSD_m <- matrix(0,l_sigmavec,l_gammavec)
PopSS_m <- matrix(0,l_sigmavec,l_gammavec)
PopSD_m <- matrix(0,l_sigmavec,l_gammavec)

res_vuln_m <- matrix(0,l_sigmavec,l_gammavec)
pop_vuln_m <- matrix(0,l_sigmavec,l_gammavec)

for (i in 1:l_sigmavec) {
  
  for (j in 1:l_gammavec) {
    
    state <- c(
      R = 0.5,            #Resource
      X = 0.5,          #Starvers
      Y = 0.5)        #Non-Starvers
    
    parameters <- c(
      alpha = 1, 
      epsilon = 1,
      sigma = sigmavec[i],
      rho = gammavec[j],
      gamma = gammavec[j],    
      mu = 0.02
    )
    t_term <- 2000
    step <- 0.1
    vuln_thresh <- 0.05
    
    
    time <- seq(0,t_term, by = step)
    
    out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
    
    ResSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),2])
    ResSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),2])
    
    XSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),3])
    XSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),3])
    YSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),4])
    YSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),4])
    
    PopSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),3] + out[((t_term/step)-1000):((t_term/step)),4])
    PopSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),3] + out[((t_term/step)-1000):((t_term/step)),4])
    
    #Vulnerability calculation
    
    res_vuln_m[i,j] <- length(which(out[,2] <= vuln_thresh))/length(out[,2])
    pop_vuln_m[i,j] <- length(which((out[,3]+out[,4]) <= vuln_thresh))/length(out[,2])
    
  }
  print(paste("i=",i))
}


PopCV_m <- PopSD_m/PopSS_m
ResCV_m <- ResSD_m/ResSS_m

pal <- wes_palette(name = "Zissou",20, type = "continuous")
M <-  PopSD_m      #(1+res_vuln_m)/(1+pop_vuln_m)
filled_contour(sigmavec, 
               gammavec, 
               M,
               levels = seq(min(M), max(M),length.out=20),col = pal,
               lwd = 0.1,xlab="p",ylab="gamma")
HopfData <- read.csv("HopfData.csv") #Import analytical solution to the Hopf Bifurcation
lines(HopfData,lwd=3,col="white")



###########################################
###########################################
###########################################



#Simulate dynamics over BOTH p and EPSILON
source("R/starv_forage_ode.R")

sigmavec <- seq(0,0.5,0.05)
l_sigmavec <- length(sigmavec)
epsilonvec <- seq(0.1,1,0.1)
l_epsilonvec <- length(epsilonvec)

ResSS_m <- matrix(0,l_sigmavec,l_epsilonvec)
ResSD_m <- matrix(0,l_sigmavec,l_epsilonvec)
XSS_m <- matrix(0,l_sigmavec,l_epsilonvec)
XSD_m <- matrix(0,l_sigmavec,l_epsilonvec)
YSS_m <- matrix(0,l_sigmavec,l_epsilonvec)
YSD_m <- matrix(0,l_sigmavec,l_epsilonvec)
PopSS_m <- matrix(0,l_sigmavec,l_epsilonvec)
PopSD_m <- matrix(0,l_sigmavec,l_epsilonvec)

res_vuln_m <- matrix(0,l_sigmavec,l_epsilonvec)
pop_vuln_m <- matrix(0,l_sigmavec,l_epsilonvec)

for (i in 1:l_sigmavec) {
  
  for (j in 1:l_epsilonvec) {
    
    state <- c(
      R = 0.5,            #Resource
      X = 0.5,          #Starvers
      Y = 0.5)        #Non-Starvers
    
    parameters <- c(
      alpha = 1, 
      epsilon = epsilonvec[j],
      sigma = sigmavec[i],
      rho = 1,
      gamma = 1,    
      mu = 0.02
    )
    t_term <- 2000
    step <- 0.1
    vuln_thresh <- 0.05
    
    
    time <- seq(0,t_term, by = step)
    
    out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
    
    ResSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),2])
    ResSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),2])
    
    XSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),3])
    XSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),3])
    YSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),4])
    YSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),4])
    
    PopSS_m[i,j] <- mean(out[((t_term/step)-1000):((t_term/step)),3] + out[((t_term/step)-1000):((t_term/step)),4])
    PopSD_m[i,j] <- sd(out[((t_term/step)-1000):((t_term/step)),3] + out[((t_term/step)-1000):((t_term/step)),4])
    
    #Vulnerability calculation
    
    res_vuln_m[i,j] <- length(which(out[,2] <= vuln_thresh))/length(out[,2])
    pop_vuln_m[i,j] <- length(which((out[,3]+out[,4]) <= vuln_thresh))/length(out[,2])
    
  }
  print(paste("i=",i))
}


PopCV_m <- PopSD_m/PopSS_m
ResCV_m <- ResSD_m/ResSS_m

pal <- wes_palette(name = "Zissou",20, type = "continuous")
M <-  PopSD_m      #(1+res_vuln_m)/(1+pop_vuln_m)
filled_contour(sigmavec, 
               epsilonvec, 
               M,
               levels = seq(min(M), max(M),length.out=20),col = pal,
               lwd = 0.1,xlab="p",ylab="gamma")
HopfData <- read.csv("HopfData.csv") #Import analytical solution to the Hopf Bifurcation
lines(HopfData,lwd=3,col="white")

