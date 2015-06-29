starvingRW <- function(L,
                       size,
                       s_max,
                       s_crit,
                       gain,
                       pr_grow,
                       pr_rep,
                       pr_mort
                       nrw,
                       srw,
                       rwloc,
                       r,
                       tmax) {
  source("R/ipbc.R")
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
  
  Rout <- list()
  Rout[[1]] <- pop_r
  Rout[[2]] <- pop_c
  
  return(Rout)
  
}