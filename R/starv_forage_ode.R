starv_forage_ode <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    
    #Rate of change
    dR <- alpha*R*(1-R) - epsilon*R*(X+Y)
    dX <- sigma*(1-R)*Y - rho*epsilon*R*X - mu*X
    dY <- gamma*epsilon*Y*R + rho*epsilon*R*X - sigma*(1-R)*Y #- mu*Y
    
    #Return rate of change
    list(c(dR,dX,dY))
  }) #end with statement
  
}