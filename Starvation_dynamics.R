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
  alpha = 0.5, 
  sigma = 0.9,
  rho = 0.2,
  lambda = 0.52,    
  mu = 0.2
)

time <- seq(0,1000, by = 0.1)

out <- ode(y = state, times = time, func = starv_forage_ode, parms = parameters)
#Plot ALL
colors <- brewer.pal(5,"Set1")
plot(out[,1],out[,2], xlab = "time", ylab = "-",type="l",lwd=2,col=colors[2],ylim=c(0,max(as.numeric(cbind(out[,2],out[,3],out[,4]))))) #Blue
lines(out[,1],out[,3],col=colors[5],type="l",lwd=2) #Orange
lines(out[,1],out[,4],col=colors[3],type="l",lwd=2) #Green

#Plot R/X+Y vs. lambdaC
with(as.list(parameters),{
  R <- out[,2]
  X <- out[,3]
  Y <- out[,4]
  C <- X + Y
  lambdaC <- lambda*Y*R - 2*mu*X - mu*Y
  ericy <- -2*sigma*Y*(1-R) + 2*rho*epsilon*R*X
  #plot((R/C),lambdaC,pch=16,cex=0.4)
  plot((R/C),ericy,pch=16,cex=0.4)
}
)

plot3d(out[,2],out[,3],out[,4],xlab="",ylab="",zlab="",type="l",col=paste(colors[2],"50",sep=""),lwd=2,xaxt="n",yaxt="n",zaxt="n",axes=FALSE)



#Simulate dynamics over BOTH sigma and lambda
source("R/starv_forage_ode.R")

sigmavec <- seq(0.0,1,0.1)
l_sigmavec <- length(sigmavec)
lambdavec <- seq(0.0,2,0.2)
l_lambdavec <- length(lambdavec)

ResSS_m <- matrix(0,l_sigmavec,l_lambdavec)
ResSD_m <- matrix(0,l_sigmavec,l_lambdavec)
XSS_m <- matrix(0,l_sigmavec,l_lambdavec)
XSD_m <- matrix(0,l_sigmavec,l_lambdavec)
YSS_m <- matrix(0,l_sigmavec,l_lambdavec)
YSD_m <- matrix(0,l_sigmavec,l_lambdavec)
PopSS_m <- matrix(0,l_sigmavec,l_lambdavec)
PopSD_m <- matrix(0,l_sigmavec,l_lambdavec)

res_vuln_m <- matrix(0,l_sigmavec,l_lambdavec)
pop_vuln_m <- matrix(0,l_sigmavec,l_lambdavec)

for (i in 1:l_sigmavec) {
  for (j in 1:l_lambdavec) {
    state <- c(
      R = 0.5,            #Resource
      X = 0.5,          #Starvers
      Y = 0.5)        #Non-Starvers
    
    parameters <- c(
      alpha = 0.5, 
      sigma = sigmavec[i],
      rho = 0.2,
      lambda = lambdavec[j],    
      mu = 0.2
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
par(mar=c(4,4,1,1))
M <-  PopSD_m      #(1+res_vuln_m)/(1+pop_vuln_m)
filled_contour(sigmavec, 
               lambdavec, 
               M,
               levels = seq(min(M), max(M),length.out=20),col = pal,
               lwd = 0.1,xlab="sigma",ylab="lambda")
mtext(expression(paste("Starvation rate, ", sigma)), side = 1, outer = F , line = 2.5)
mtext(expression(paste("Consumer growth rate, ", lambda)), side = 2, outer = F , line = 2.5)
HopfData <- read.csv("HopfDataFR.csv",header=FALSE) #Import analytical solution to the Hopf Bifurcation
lines(HopfData,col="white",lwd=3)
#SData <- cbind(seq(0.1,1,0.01),seq(0.1,1,0.01))
#lines(SData,col=colors[2],lwd=3)



