library(deSolve)
#library(rgl)
library(RColorBrewer)
library(wesanderson)
source("R/filled_contour.r")
source("R/smooth_pal.R")

#Starving forager ODE
source("R/comp_forage_ode.R")

state <- c(
  R = 0.5,            
  X1 = 0.5,          
  Y1 = 0.5,
  X2 = 0.5,
  Y2 = 0.5)        

parameters <- c(
  alpha = 1, 
  epsilon = 0.5,
  sigma1 = 0.2,
  sigma2 = 0.1,
  rho = 0.2,
  gamma = 0.8,    
  mu = 0.01
)

time <- seq(0,1000, by = 0.1)

out <- ode(y = state, times = time, func = comp_forage_ode, parms = parameters)
#Plot ALL
colors <- brewer.pal(5,"Set1")
plot(out[,1],out[,2], xlab = "time", ylab = "-",type="l",lwd=2,col=colors[2],
     ylim=c(0,max(as.numeric(cbind(out[,2],out[,3],out[,4],out[,5],out[,6]))))) #Blue
lines(out[,1],out[,3],col=colors[5],type="l",lwd=2) #Orange
lines(out[,1],out[,4],col=colors[3],type="l",lwd=2) #Green
lines(out[,1],out[,5],col=colors[5],type="l",lwd=2,lty=2) #Orange
lines(out[,1],out[,6],col=colors[3],type="l",lwd=2,lty=2) #Green


#Who wins sigma1 vs. sigma2
source("R/comp_forage_ode.R")
sigma1_vec <- seq(0,1,0.1)
sigma2_vec <- seq(0,1,0.1)
rho1_vec <- sigma1_vec
rho2_vec <- sigma2_vec
l_sigma1 <- length(sigma1_vec)
l_sigma2 <- length(sigma2_vec)
pr_m <- matrix(0,l_sigma1,l_sigma2)
pr_c <- matrix(0,(l_sigma1*l_sigma2),3)
tic <- 0
for (i in 1:l_sigma1) {
  for (j in 1:l_sigma2) {
    tic <- tic + 1
    state <- c(
      R = 0.5,            
      X1 = 0.5,          
      Y1 = 0.5,
      X2 = 0.5,
      Y2 = 0.5)        
    parameters <- c(
      alpha = 1, 
      epsilon = 0.5,
      sigma1 = sigma1_vec[i],
      sigma2 = sigma2_vec[j],
      rho1 = rho1_vec[i],
      rho2 = rho2_vec[j],
      gamma = 0.8,    
      mu = 0.01
    )
    time <- seq(0,200, by = 0.1)
    out <- ode(y = state, times = time, func = comp_forage_ode, parms = parameters)
    pop1 <- out[,3] + out[,4]
    pop2 <- out[,5] + out[,6]
    pop_ratio <- pop1/pop2
    pr_m[i,j] <- median(pop_ratio[(length(pop_ratio)-100):(length(pop_ratio))])
    pr_c[tic,] <- c(i,j,pr_m[i,j])
  }
}

bw <- pr_c[,3] > 1

plot(pr_c[,1],pr_c[,2],pch=16,col=bw)
points(pr_c[,1],pr_c[,2])

