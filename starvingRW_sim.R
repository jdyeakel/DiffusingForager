library(RColorBrewer)
library(Rcpp)
library(animation)
#nearest neighbor with boundary conditions function
source("R/ipbc.R")

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

############
#Parameters
############

#Probability of resource growth | presence of nearest neighbors
pr_grow <- 0.5
#Probability of consumer reproduction
pr_rep <- 0.5
#Probability of consumer mortality | they are starving
pr_mort <- 1


#initial number of random walkers
nrw <- 10
#Random walker energetic state
srw <- rep(s_max,nrw)
#Location of RWs
rwloc <- sample(seq(1:size),nrw,replace=TRUE)
#Resource state at each location
r <- rep(1,size)

#Maximum time
tmax <- 100

#Vectors for consumer and resource population sizes over time
pop_c <- numeric(tmax)
pop_r <- numeric(tmax)
r_frame <- list()
for (t in 1:tmax) {
  
  #Across each individual in the system...
  ind_check <- TRUE
  num <- length(srw)
  pop_c[t] <- num
  pop_r[t] <- sum(r)
  r_frame[[t]] <- r
  #rwloc_frame[[t]] <- rwloc
  i <- 0
  
  #Across the number of individuals (which will fluctuate)
  while (ind_check == TRUE) {
    
    #individual index
    i <- i + 1
    
    #Move to a new location
    new_loc <- sample(ipbc(rwloc[i],L),1)
    #new_loc <- sample(seq(1,size),1)
    #consume resource if it is there
    new_s <- min(srw[i] - 1 + r[new_loc]*gain,s_max)
    #deplete the resource
    r[new_loc] <- 0
    
    #Impliment mortality and update
    #Random mortality
    if (new_s <= 0) {
      #Draw random value
      rdraw <- runif(1)
      if (rdraw < pr_mort) {
        #Spaghetti method of removal
        srw[i] <- srw[num]; srw <- srw[1:(num-1)]
        rwloc[i] <- rwloc[num]; rwloc <- rwloc[1:(num-1)]
        i <- i - 1 #reanalyzes the 'new member' of the slot
      } else {
        #Individual survives
        #Record the new location and new state
        rwloc[i] <- new_loc
        srw[i] <- new_s
      }
    } else {
      #Individual survives:
      #Record the new location and new state
      rwloc[i] <- new_loc
      srw[i] <- new_s
    }
    
    #Recalculate num (number of individuals)... will be shorter if there is mortality
    num <- length(srw)
    if (i == num) {
      ind_check <- FALSE #The while loop is stopped if you get to the end of the list
    }
    
  } #End while loop
  
  #4) RW reproduction
  for (i in 1:num) {
    #Reproduction can only happen if s > s_crit
    if (srw[i] > s_crit) {
      #Draw random value
      rdraw <- runif(1)
      if (rdraw < pr_rep) {
        #Add a new individual to the end of the vector
        srw_new <- srw[i]
        rwloc_new <- rwloc[i]
        #rwloc_new <- sample(seq(1:size),1)
        srw <- c(srw,srw_new)
        rwloc <- c(rwloc,rwloc_new)
        #Recalculate the number of individuals... will be greater if there is reproduction
        num <- length(srw)
      }
    }
  }
  
  #Resource growth
  for (j in 1:length(r)) {
    if (r[j] == 0) {
      #Determine the number of nearest neighbor resources
      r_nn <- sum(r[ipbc(j,L)])
      #r_nn <- r[sample(seq(1:size),1)]
      #Growth only occurs if there is at least one nearest neighbor
      if (r_nn >= 1) {
        r_draw <- runif(1)
        if (r_draw < pr_grow) {
          r[j] <- 1
        }
      }
    }
  }
 #Stop the simulation if the number of individuals is equal to one... never gets to zero - not sure why.
  if (num == 1) {
    break
  }
  
  
}

sourceCpp("src/starvingRW.cpp")
cout <- starvingRW(L, s_max, s_crit, gain, tmax, pr_grow, pr_rep, pr_mort, srw, rwloc-1, r)
pop_r <- cout[[1]]
pop_c <- cout[[2]]
r_frame <- cout[[3]]
srwm <- cout[[4]]

pal <- brewer.pal(3,"Set1")
plot(pop_c,ylim=c(0,(L+2)^2),type="l",col=pal[1],lwd=2)
points(pop_r,type="l",col=pal[2],lwd=2)

plot(pop_r,pop_c,pch=16,col=pal[2])

#http://programmingr.com/content/animations-r/
ani.options(interval=.1)
saveGIF({
  for (i in 1:300) {
    image(matrix(r_frame[[i]],(L+2),(L+2)),col=c("white","black"))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/animations/resource.gif")


ani.options(interval=.1)
saveGIF({
  for (i in 1:tmax) {
    rl <- r_frame[i,]
    image(matrix(rl,(L+2),(L+2)),col=c("white","black"))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/animations/resource.gif")




