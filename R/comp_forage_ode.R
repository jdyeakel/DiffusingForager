comp_forage_ode <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    
    #Rate of change
    dR <- alpha*R*(1-R) - epsilon*R*(X1+Y1) - epsilon*R*(X2+Y2)
    dX1 <- sigma1*(1-R)*Y1 - rho1*epsilon*R*X1 - mu*X1
    dY1 <- gamma1*epsilon*Y1*R + rho1*epsilon*R*X1 - sigma1*(1-R)*Y1 #- mu*Y1
    dX2 <- sigma2*(1-R)*Y2 - rho2*epsilon*R*X2 - mu*X2
    dY2 <- gamma2*epsilon*Y2*R + rho2*epsilon*R*X2 - sigma2*(1-R)*Y2 #- mu*Y2
    
    #Return rate of change
    list(c(dR,dX1,dY1,dX2,dY2))
  }) #end with statement
  
}