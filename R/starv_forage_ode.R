starv_forage_ode <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    
    #Rate of change
    dR <- (alpha+rnorm(1,0,epsilon))*R*(1-R) - R*(rho*H+beta*Full)
    dH <- sigma*(1-R)*Full - rho*R*H - mu*H
    dFull <- lambda*Full + rho*R*H - sigma*(1-R)*Full
    
    #Return rate of change
    list(c(dR,dH,dFull))
  }) #end with statement
  
}