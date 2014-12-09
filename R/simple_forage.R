simple_forage <- function(dim,xmax,dist_max,rep) {
  
  tmort <- numeric(rep)
  path_list <- list()
  x_list <- list()
  
  for (r in 1:rep) {
    
    origin <- numeric(dim)
    res_array <- array(data = 1, dim = rep(dist_max,dim))
    
    path <- numeric(dim)
    x <- numeric()
    
    x[1] <- xmax
    t = 1
    while (x[t] != 0) {
      t = t+1
      #Define the next step of the consumer's path
      step <- sample(c(-1,0,1),dim,replace=TRUE)
      
      if (t == 2) {
        newpath <- path + step
      } else {
        newpath <- path[t-1,] + step
      }
      
      
      #Boundary conditions
      toolow <- which(newpath < 1)
      if (length(toolow) > 0) {
        newpath[toolow] <- 1
      }
      
      #accept proposed step
      path <- rbind(path,newpath)
      
      x[t] <- max(x[t-1] - 1, 0)
      
      if (res_array[newpath] == 1) {
        x[t] <- xmax
      }
      
      #Eliminate resource
      res_array[newpath] <- 0
      
    } #end t
    
    tmort[r] <- t
    path_list[[r]] <- path
    x_list[[r]] <- x
    #out <- list()
    #out[[1]] <- path
    #out[[2]] <- x
    
    #return(out)
  } #end rep
  
  out <- list()
  out[[1]] <- tmort
  out[[2]] <- path_list
  out[[3]] <- x_list
  
  return(out)
}