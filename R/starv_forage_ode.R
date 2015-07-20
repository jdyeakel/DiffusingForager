starv_forage_ode <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    
    #Rate of change
    dR <- alpha*R*(1-R) - R*(X+Y)
    dX <- sigma*(1-R)*Y - rho*R*X - mu*X
    dY <- lambda*Y + rho*R*X - sigma*(1-R)*Y
    
    #Return rate of change
    list(c(dR,dX,dY))
  }) #end with statement
  
}