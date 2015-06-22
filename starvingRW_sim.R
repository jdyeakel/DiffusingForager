

#Cell lattice
#How many cells? (total)
L <- 10
size <- (L+2)^2
nloc <- seq(1,size,1)
#Resource state at each location
r <- rep(1,size)

ipbc <- function(x,L) {
  check <- 0
  if (x <= L+2) {new_x <- x + L*(L+2); check <- 1}
  if (x >= ((L+2)^2 - (L+2))) {new_x <- x - L*(L+2); check <- 1}
  if (x%%(L+2) == 1) {new_x <- x + (L); check <- 1}
  if (x%%(L+2) == 0) {new_x <- x - (L); check <- 1}
  if(check == 0) {new_x <- x}
  nn <- c(new_x+1,new_x-1,new_x+(L+2),new_x-(L+2))
  return(nn)
}

#Maximum energetic state
s_max <- 50
s_crit <- 25
gain <- s_max

#initial number of random walkers
nrw <- 10
#Random walker energetic state
srw <- rep(s_max,nrw)

rwloc <- sample(nloc,x,replace=TRUE)

tmax <- 1000

#Should all be written in CPP
for (t in 1:tmax) {
  
  #6) RW Move
  rwloc <- function(rwloc) {}
  
  # Internal Dynamics
  
  #1) update s vector and r vector
  srw <- srw - 1 + r[rwloc]*gain # Minus 1 for prior movement, and + resource for consumption
  #Boundary conditions... can't be greater than s_max
  srw <- sapply(srw,function(x){min(s_max,x)})
  
  #2) update r vector
  r[rwloc] <- r[rwloc] - 1
  
  #3) implement starvation and mortality
  dead <- which(srw <= 0)
  srw <- srw[-dead]
  
  
  
  #4) RW reproduction
  full <- which(srw > s_crit) #Alternatively, write reproduction as a probability as a function of full
  srw_new <- srw[full] #begin with same S as parents
  srw <- c(srw,srw_new)
  rwloc_new <- rwloc[full] #begin at same cell location as parents
  rwloc <- c(rwloc,rwloc_new)
  
  
  #5) Resource growth ~ from NEIGHBORS
  #The number of nearest neighbors for each resource patch
  r_nn <- function(r) {}
  # initiate regrowth
  r_grow <- which(r_nn > 0) #Alternatively, write colonization as a probability as a funciton of r_nn
  r[r_grow] <- 1

  #mod = %%

  
}