

#Cell lattice
#How many cells?
n <- 1000
nloc <- seq(1,n,1)
#Resource state at each location
r <- rep(1,n)
#Resource delay clock
r_tic <- rep(0,n)
#delay time
t_delay <- 5

#Maximum energetic state
s_max <- 50
s_crit <- 25
#initial number of random walkers
nrw <- 10
#Random walker energetic state
srw <- rep(s_max,nrw)

#Set record keeping slots
rwloc <- sample(nloc,x,replace=TRUE)

tmax <- 1000

for (t in 1:tmax) {
  
  # Internal Dynamics
  #1) update s vector and r vector
  srw <- srw - 1 + r[rwloc]
  #2) update r vector
  r[rwloc] <- r[rwloc] - 1
  
  #3) implement starvation and mortality
  dead <- which(srw == 0)
  srw <- srw[-dead]
  
  
  
  # RW reproduction
  full <- which(srw > s_crit)
  srw_new <- srw[full]
  srw <- c(srw,srw_new)
  rwloc_new <- rwloc[full]
  rwloc <- c(rwloc,rwloc_new)
  
  
  # Resource growth and colonization
  r_tic <- r_tic + 1
  # initiate regrowth
  r_grow <- which(r_tic == t_delay)
  r[r_grow] <- 1
  #Reset clock
  r_tic[r_grow] <- 0
  
  # RW Move
  rwloc <- FUNCTION
  
}