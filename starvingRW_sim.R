library(RColorBrewer)
library(Rcpp)
library(animation)
#nearest neighbor with boundary conditions function
source("R/ipbc.R")
sourceCpp("src/starvingRW_pr.cpp")

#Cell lattice
#How many cells? (total)
L <- 50
size <- (L+2)^2
#nloc <- seq(1,size,1)

#
#RW new location rand
#Offspring location rand
#Resource nn rand


#Maximum energetic state
s_max <- 1
s_crit <- 0
gain <- s_max

#initial number of random walkers
nrw <- floor(0.1*size)
#Random walker energetic state
srw <- rep(s_max,nrw)
#Location of RWs
rwloc <- sample(seq(1:size),nrw,replace=TRUE)
#Resource state at each location
r <- rep(1,size)

#If p = 0, we have the lattice dynamic
#If p = 1, we have the mean field dynamic
p <- 1

#Maximum time
tmax <- 1000

############
#Parameters
############

#Probability of resource growth | presence of nearest neighbors
alpha <- 0.5
#Rate of starvation
sigma <- 1
#Rate of recovery
rho <- 0.2
#Probability of consumer reproduction
lambda <- 0.1
#Probability of consumer mortality | they are starving
mu <- 0.2


cout <- starvingRW_pr(L, s_max, s_crit, gain, tmax, alpha, sigma, rho, lambda, mu, srw, rwloc-1, r, p)
pop_r <- cout[[1]]
pop_c <- cout[[2]]
r_frame <- cout[[3]]
locm <- cout[[4]]
srwm <- cout[[5]]
propr <- pop_r/size #(pop_r + pop_c)
propc <- pop_c/size #(pop_r + pop_c)

pal <- brewer.pal(3,"Set1")
plot(pop_c,ylim=c(0,max(c(pop_c,pop_r))),type="l",col=pal[1],lwd=2)
points(pop_r,type="l",col=pal[2],lwd=2)


plot(propc,ylim=c(0,1),type="l",col=pal[1],lwd=2)
points(propr,type="l",col=pal[2],lwd=2)


plot(pop_r,pop_c,pch=".",col=pal[2])
plot(diff(pop_r),diff(pop_c),pch=".",col=pal[2])


ani.options(interval=.025)
saveGIF({
  for (i in 1:tmax) {
    rl <- r_frame[i,]
    image(matrix(rl,(L+2),(L+2)),col=c("white","black"))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/animations/resourceL100_rand.gif")



#Fixed point analysis


#Probability of resource growth | presence of nearest neighbors
alpha <- 0.5
#Rate of starvation
#sigma <- 1
#Rate of recovery
rho <- 0.2
#Probability of consumer reproduction
lambda <- 0.1
#Probability of consumer mortality | they are starving
mu <- 0.2
tmax <- 1000
sigmaseq <- seq(0.1,1,0.1)
l_sigmaseq <- length(sigmaseq)
meanr <- numeric(l_sigmaseq)
meanc <- numeric(l_sigmaseq)
meanpropr <- numeric(l_sigmaseq)
meanpropc <- numeric(l_sigmaseq)
tic <- 0
for (sigma in sigmaseq) {
  tic <- tic + 1
  print(tic)
  cout <- starvingRW_pr(L, s_max, s_crit, gain, tmax, alpha, sigma, rho, lambda, mu, srw, rwloc-1, r, p)
  pop_r <- cout[[1]]
  pop_c <- cout[[2]]
  propr <- pop_r/(pop_r + pop_c)
  propc <- pop_c/(pop_r + pop_c)
  meanr[tic] <- mean(pop_r[floor(tmax/2):1000])
  meanc[tic] <- mean(pop_c[floor(tmax/2):1000])
  meanpropr[tic] <- mean(propr[floor(tmax/2):1000])
  meanpropc[tic] <- mean(propc[floor(tmax/2):1000])
}
plot(meanr,type="l",col=pal[2],ylim=c(0,size))
lines(meanc,col=pal[1])

plot(sigmaseq,meanpropr,type="l",col=pal[2],ylim=c(0,1),xlim=c(0,1))
lines(sigmaseq,meanpropc,col=pal[1])






