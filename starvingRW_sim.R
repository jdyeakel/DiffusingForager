library(RColorBrewer)
library(wesanderson)
library(Rcpp)
library(animation)
#nearest neighbor with boundary conditions function
source("R/ipbc.R")
sourceCpp("src/starvingforager_eventNM.cpp")

#Initiate starting conditions
L <- 50
size <- (L-2)^2
t_term <- 10

#Parameters
alpha <- 1
K <- size
sigma <- 0.5
rho <- 0.2
lambda <- 0.3
mu <- 0.2
#D <- 1


# ind_vec <- sample(c(0,1,2),size,replace=T)
# loc_vec <- sample(seq(0,size-1),size,replace=T)

cons_init <- 1000
ind_vec <- c(numeric(size),sample(c(1,2),cons_init,replace=T))
loc_vec <- c(sample(seq(0,size-1),size,replace=T),sample(seq(0,size-1),cons_init,replace=T))

Rout <- starvingforager_eventNM(
  L,
  t_term,
  alpha,
  K,
  sigma,
  rho,
  lambda,
  mu,
  #D,
  ind_vec,
  loc_vec)

state <- Rout[[1]]
loc <- Rout[[2]]
time <- Rout[[3]]

Rden <- unlist(lapply(state,function(x){length(which(x == 0))})) #/size
Sden <- unlist(lapply(state,function(x){length(which(x == 1))})) #/size
Fden <- unlist(lapply(state,function(x){length(which(x == 2))})) #/size

pal <- brewer.pal(5,"Set1")
plot(time,Rden,type="l",ylim=c(0,max(Rden,Sden,Fden)),col=pal[2],lwd=3)
lines(time,Sden,col=pal[5],lwd=3)
lines(time,Fden,col=pal[3],lwd=3)














#Cell lattice
#How many cells? (total)
L <- 200
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
tmax <- 5000

############
#Parameters
############

#Probability of resource growth | presence of nearest neighbors
alpha <- 0.5
#Rate of starvation
sigma <- 0.12
#Rate of recovery
rho <- 0.2
#Probability of consumer reproduction
lambda <- 0.1
#Probability of consumer mortality | they are starving
mu <- 0.2


cout <- starvingRW_pr(L, s_max, s_crit, gain, tmax, alpha, sigma, rho, lambda, mu, srw, rwloc-1, r, p)
pop_r <- cout[[1]]
pop_c <- cout[[2]]
pop_s <- cout[[3]]
pop_f <- cout[[4]]
r_frame <- cout[[5]]
locm <- cout[[6]]
srwm <- cout[[7]]
propr <- pop_r/size #(pop_r + pop_c)
propc <- pop_c/size #(pop_r + pop_c)

pal <- brewer.pal(5,"Set1")
plot(pop_r,ylim=c(0,max(c(pop_c,pop_r))),type="l",col=pal[2],lwd=2)
lines(pop_s,type="l",col=pal[5],lwd=2)
lines(pop_f,type="l",col=pal[3],lwd=2)


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
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/animations/resourceL50_pulse.gif")



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
sigmaseq <- seq(0.12,1,0.1)
l_sigmaseq <- length(sigmaseq)
meanr <- numeric(l_sigmaseq)
meanc <- numeric(l_sigmaseq)
means <- numeric(l_sigmaseq)
meanf <- numeric(l_sigmaseq)
meanpropr <- numeric(l_sigmaseq)
meanpropc <- numeric(l_sigmaseq)
meanprops <- numeric(l_sigmaseq)
meanpropf <- numeric(l_sigmaseq)
tic <- 0
for (sigma in sigmaseq) {
  tic <- tic + 1
  print(tic)
  cout <- starvingRW_pr(L, s_max, s_crit, gain, tmax, alpha, sigma, rho, lambda, mu, srw, rwloc-1, r, p)
  pop_r <- cout[[1]]
  pop_c <- cout[[2]]
  pop_s <- cout[[3]]
  pop_f <- cout[[4]]
  propr <- pop_r/(pop_r + pop_c)
  propc <- pop_c/(pop_r + pop_c)
  props <- pop_s/(pop_r + pop_c)
  propf <- pop_f/(pop_r + pop_c)
  meanr[tic] <- mean(pop_r[floor(tmax/2):1000])
  meanc[tic] <- mean(pop_c[floor(tmax/2):1000])
  means[tic] <- mean(pop_s[floor(tmax/2):1000])
  meanf[tic] <- mean(pop_f[floor(tmax/2):1000])
  meanpropr[tic] <- mean(propr[floor(tmax/2):1000])
  meanpropc[tic] <- mean(propc[floor(tmax/2):1000])
  meanprops[tic] <- mean(props[floor(tmax/2):1000])
  meanpropf[tic] <- mean(propf[floor(tmax/2):1000])
  print(paste("sigma = ",sigma))
}
# plot(meanr,type="l",col=pal[2],ylim=c(0,size))
# lines(means, type="l",col=pal[5])
# lines(meanf, type="l",col=pal[3])
# lines(meanc,col=pal[1])

plot(sigmaseq,meanpropr,type="l",col=pal[2],
     ylim=c(0,1),xlim=c(0,1),lwd=3,ylab="Prop. Abundance",xlab=expression(paste("Starvation rate ",sigma)))
lines(sigmaseq,meanprops, type="l",col=pal[5],lwd=3)
lines(sigmaseq,meanpropf, type="l",col=pal[3],lwd=3)
lines(sigmaseq,meanpropc,col=pal[1])


#Sigma vs. lambda
#System Size
L <- 50
size <- (L+2)^2
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
#Probability of resource growth | presence of nearest neighbors
alpha <- 0.5
#Rate of recovery
rho <- 0.2
#Probability of consumer mortality | they are starving
mu <- 0.2
#Spatial (p=0) vs. Mean Field (p=1)
p <- 0

tmax <- 500
res_burn_l <- list()
c_burn_l <- list()
lambdaseq <- seq(0.1,1,0.05)
sigmaseq <- seq(0.1,1,0.05)
r_sd <- matrix(0,length(sigmaseq),length(lambdaseq))
c_sd <- matrix(0,length(sigmaseq),length(lambdaseq))
toc <- 0
pop_r_traj <- list()
pop_c_traj <- list()
for (i in 1:length(lambdaseq)) {
  pop_r_traj[[i]] <- list()
  pop_c_traj[[i]] <- list()
}

for (lambda in lambdaseq) {
  toc <- toc + 1
  #   res_burn_s <- list()
  #   c_burn_s <- list()
  tic <- 0
  for (sigma in sigmaseq) {
    tic <- tic + 1
    if (lambda < sigma) {
      #Initialize the state for error-reading
      cout <- "resource extinction"
      #If there is an out of bounds error in the cpp code, it will be rerun
      extinct_tic <- 0
      print(paste("lambda=",lambda,"   sigma=",sigma))
      while((cout[[1]] == "consumer extinction") || (cout[[1]] == "resource extinction")) { 
        cout <- starvingRW_pr(L, s_max, s_crit, gain, tmax, alpha, sigma, rho, lambda, mu, srw, rwloc-1, r, p)
        pop_r <- cout[[1]]
        pop_c <- cout[[2]]
        r_sd[tic,toc] <- sd(pop_r[floor(tmax/2):tmax])
        c_sd[tic,toc] <- sd(pop_c[floor(tmax/2):tmax])
        #Save trajectories!
        pop_r_traj[[toc]][[tic]] <- pop_r
        pop_c_traj[[toc]][[tic]] <- pop_c
        extinct_tic <- extinct_tic + 1
        print(paste("extinctions = ",extinct_tic))
        if ((extinct_tic == 20) && (cout[[1]] == "resource extinction")) {
          r_sd[tic,toc] <- -1
          c_sd[tic,toc] <- -1
          pop_r_traj[[toc]][[tic]] <- -1
          pop_c_traj[[toc]][[tic]] <- -1
          break
        }
        if ((extinct_tic == 20) && (cout[[1]] == "consumer extinction")) {
          r_sd[tic,toc] <- -2
          c_sd[tic,toc] <- -2
          pop_r_traj[[toc]][[tic]] <- -2
          pop_c_traj[[toc]][[tic]] <- -2
          break
        }
        if ((extinct_tic == 20) && (cout[[1]] == "consumer overflow")) {
          r_sd[tic,toc] <- -3
          c_sd[tic,toc] <- -3
          pop_r_traj[[toc]][[tic]] <- -3
          pop_c_traj[[toc]][[tic]] <- -3
          break
        }
      }
      print(paste("toc= ",toc,"; tic= ",tic))
    } 
  }
  #   res_burn_l[[toc]] <- res_burn_s
  #   c_burn_l[[toc]] <- res_burn_s
}
data <- list(); data[[1]] <- pop_r_traj;data[[2]] <- pop_c_traj
#save(data,file="/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/lambda_sigma_traj.RData")

library(plotrix)
mfreq <- matrix(0,length(lambdaseq),length(sigmaseq))
pal <- matrix(0,length(lambdaseq),length(sigmaseq))
colornames <- smooth_pal(brewer.pal(11,"Spectral"),n=4)
colorvalues <- seq(0,0.04,0.001)
for (i in 1:length(lambdaseq)) {
  for (j in 1:length(sigmaseq)) {
    if (i < j) {
      r_traj <- pop_r_traj[[i]][[j]]
      #c_traj <- pop_c_traj[[i]][[j]]
      if (length(r_traj) > 5) {
        r_traj <- r_traj[250:500]
        if (max(r_traj) == min(r_traj)) {
          mfreq[i,j] <- 0
          colorpos <- which.min(abs(colorvalues - mfreq[i,j])) 
          pal[i,j] <- colornames[colorpos]
        } else {
          sv <- spectrum(r_traj);
          mfreq[i,j] <- sv$freq[which(sv$spec == max(sv$spec))]
          colorpos <- which.min(abs(colorvalues - mfreq[i,j])) 
          pal[i,j] <- colornames[colorpos]
        }
      } else {
        #The case of resource extinction
        if (r_traj == -1) {
          mfreq[i,j] <- -1
          pal[i,j] <- "black"
        }
        #The case of consumer extinction
        if (r_traj == -2) {
          mfreq[i,j] <- -2
          pal[i,j] <- "gray"
        }
        #The case of consumer overflow (~resource extinction)
        if ((r_traj == "consumer overflow") || (r_traj == -3)) {
          mfreq[i,j] <- -3
          pal[i,j] <- "black"
        }
      }
    } else {pal[i,j] <- "white"}
  } 
}
color2D.matplot(mfreq, border="white", axes=FALSE, cellcolors=pal,yrev=FALSE,
                xlab=expression(paste("Starvation rate  ",sigma)),
                ylab=expression(paste("Consumer growth rate  ",lambda)))
lines(seq(0,length(lambdaseq)),seq(0,length(lambdaseq)))
axis(1, at = seq(1,19,2), labels = seq(0.1,1,0.1))
axis(2, at = seq(1,19,2), labels = seq(0.1,1,0.1))

pal <- wes_palette(name = "Zissou",20, type = "continuous")
par(mar=c(4,4,1,1))
M <-  c_sd      #(1+res_vuln_m)/(1+pop_vuln_m)
filled_contour(sigmaseq, 
               lambdaseq, 
               M,
               levels = seq(min(M), max(M),length.out=20),col = pal,
               lwd = 0.1,xlab="sigma",ylab="lambda")
mtext(expression(paste("Starvation rate, ", sigma)), side = 1, outer = F , line = 2.5)
mtext(expression(paste("Consumer growth rate, ", lambda)), side = 2, outer = F , line = 2.5)
HopfData <- read.csv("HopfData.csv") #Import analytical solution to the Hopf Bifurcation
SData <- cbind(seq(0.1,1,0.01),seq(0.1,1,0.01))
lines(SData,col=colors[2],lwd=3)
lines(HopfData,col=colors[5],lwd=3)



