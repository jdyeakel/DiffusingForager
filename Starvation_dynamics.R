#Starvation dynamics

library(deSolve)
library(rgl)
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
  alpha = 1,
  k = 2,
  p = 0.8,       
  mu_x = 0.02,     
  mu_y = 0.002,
  lambda = 1.7
)

time <- seq(0,5000, by = 0.1)

out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
#Plot ALL
colors <- brewer.pal(5,"Set1")
plot(out[,1],out[,2], xlab = "time", ylab = "-",type="l",lwd=2,col=colors[2],ylim=c(0,2)) #Blue
lines(out[,1],out[,3],col=colors[5],type="l",lwd=2) #Orange
lines(out[,1],out[,4],col=colors[3],type="l",lwd=2) #Green

#Vulnerability calculation
vuln_thresh <- 0.1
res_vuln <- length(which(out[,2] <= 0.1))/length(out[,2])

#2D and 3D plots
plot(out[,2],out[,3],type="l",xlab="Resource",ylab="Foragers")
plot3d(out[,2],out[,3],out[,4],xlab="",ylab="",zlab="",type="l",col=paste(colors[2],"50",sep=""),lwd=2,xaxt="n",yaxt="n",zaxt="n")

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
    mu_x = 0.2,     
    mu_y = 0.1,
    lambda = 5
  )
  
  time <- seq(0,500, by = 0.1)
  
  out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
  
  ResSS[i] <- median(out[4000:5000,2])
  ResSD[i] <- sd(out[4000:5000,2])
  
  PopSS[i] <- median(out[4000:5000,3] + out[4000:5000,4])
  PopSD[i] <- sd(out[4000:5000,3] + out[4000:5000,4])
  
}

plot(pvec,ResSS,ylim=c(0,1),type="b")
plot(pvec,ResSD,ylim=c(0,1),type="b")
plot(pvec,ResSD/ResSS,ylim=c(0,max(ResSD/ResSS)),type="b")

plot(pvec,PopSS,ylim=c(0,1),type="b")
plot(pvec,PopSD,ylim=c(0,1),type="b")
plot(pvec,PopSD/PopSS,ylim=c(0,max(PopSD/PopSS)),type="b")

colors <- brewer.pal(5,"Set1")
plot(pvec,PopSD/PopSS,ylim=c(0,max(PopSD/PopSS,ResSD/ResSS)),type="b",pch=16,col=colors[1])
points(pvec,ResSD/ResSS,ylim=c(0,max(ResSD/ResSS)),type="b",pch=16,col=colors[2])









#Simulate dynamics over BOTH p and lambda
source("R/starv_forage_ode.R")

pvec <- seq(0,0.5,0.05)
l_pvec <- length(pvec)
lambdavec <- seq(0.1,2,0.1)
l_lambdavec <- length(lambdavec)

ResSS_m <- matrix(0,l_pvec,l_lambdavec)
ResSD_m <- matrix(0,l_pvec,l_lambdavec)
XSS_m <- matrix(0,l_pvec,l_lambdavec)
XSD_m <- matrix(0,l_pvec,l_lambdavec)
YSS_m <- matrix(0,l_pvec,l_lambdavec)
YSD_m <- matrix(0,l_pvec,l_lambdavec)
PopSS_m <- matrix(0,l_pvec,l_lambdavec)
PopSD_m <- matrix(0,l_pvec,l_lambdavec)

res_vuln_m <- matrix(0,l_pvec,l_lambdavec)
pop_vuln_m <- matrix(0,l_pvec,l_lambdavec)

for (i in 1:l_pvec) {
  
  for (j in 1:l_lambdavec) {
    
    state <- c(
      R = 0.5,            #Resource
      X = 0.5,          #Starvers
      Y = 0.5)        #Non-Starvers
    
    parameters <- c(
      alpha = 2,
      k = 1,
      p = pvec[i],       
      mu_x = 0.02,     
      mu_y = 0.002,
      lambda = lambdavec[j]
    )
    t_term <- 1000
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

pal <- wes.palette(name = "Zissou", type = "continuous")
M <-  ResSD_m      #(1+res_vuln_m)/(1+pop_vuln_m)
filled_contour(pvec, 
               lambdavec, 
               M,
               levels = seq(min(M), max(M),length.out=20),col = pal(20),
               lwd = 0.1,xlab="p",ylab="lambda")
HopfData <- read.csv("HopfData.csv") #Import analytical solution to the Hopf Bifurcation
lines(HopfData,lwd=3,col="white")



