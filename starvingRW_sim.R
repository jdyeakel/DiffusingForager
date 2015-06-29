library(RColorBrewer)
library(Rcpp)
library(animation)
#nearest neighbor with boundary conditions function
source("R/ipbc.R")

#Cell lattice
#How many cells? (total)
L <- 100
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

############
#Parameters
############

#Probability of resource growth | presence of nearest neighbors
pr_grow <- 0.01
#Probability of consumer reproduction
pr_rep <- 0.5
#Probability of consumer mortality | they are starving
pr_mort <- 0.25


#initial number of random walkers
nrw <- 10
#Random walker energetic state
srw <- rep(s_max,nrw)
#Location of RWs
rwloc <- sample(seq(1:size),nrw,replace=TRUE)
#Resource state at each location
r <- rep(1,size)

#Maximum time
tmax <- 5000

sourceCpp("src/starvingRW.cpp")
cout <- starvingRW(L, s_max, s_crit, gain, tmax, pr_grow, pr_rep, pr_mort, srw, rwloc-1, r)
pop_r <- cout[[1]]
pop_c <- cout[[2]]
r_frame <- cout[[3]]
locm <- cout[[4]]
srwm <- cout[[5]]

pal <- brewer.pal(3,"Set1")
plot(pop_c,ylim=c(0,max(c(pop_c,pop_r))),type="l",col=pal[1],lwd=2)
points(pop_r,type="l",col=pal[2],lwd=2)

plot(pop_r,pop_c,pch=".",col=pal[2])
plot(diff(pop_r),diff(pop_c),pch=".",col=pal[2])


ani.options(interval=.025)
saveGIF({
  for (i in 1:tmax) {
    rl <- r_frame[i,]
    image(matrix(rl,(L+2),(L+2)),col=c("white","black"))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/animations/resourceL200.gif")




