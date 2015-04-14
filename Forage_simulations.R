source("R/simple_forage.R")

dim <- 1
xmax <- 100
dist_max <- 1000
rep <- 100
rout <- simple_forage(dim=dim,xmax=xmax,dist_max=dist_max,rep=rep)

#Plot spatial trajectories
plot(rout[[2]][[1]],type="l",col="lightgray",xlim=c(0,max(rout[[1]])),ylim=c(0,max(unlist(rout[[2]]))))
for (i in 2:rep) {
  lines(rout[[2]][[i]],col="lightgray")
}

#Plot energetic trajectories
plot(rout[[3]][[1]],type="l",col="lightgray",xlim=c(0,max(rout[[1]])),ylim=c(0,max(unlist(rout[[3]]))))
for (i in 2:1000) {
  lines(rout[[3]][[i]],col="lightgray")
}

#Across dimensions
max.dim <- 4
m.mort <- numeric(max.dim)
for (i in 1:max.dim) {
  print(paste("i=",i,sep=""))
  rout <- simple_forage(dim=i,xmax=10,dist_max=100,rep=100)
  m.mort[i] <- mean(rout[[1]])
}