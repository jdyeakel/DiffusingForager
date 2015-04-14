starv_forage_ode <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    
    #Rate of change
    dR <- alpha*R*(1-R) - R*(X+Y)
    dX <- p*(1-R)*Y - R*X - mu_x*X
    dY <- lambda*Y*R + R*X - p*(1-R)*Y - mu_y*Y
    
    #Return rate of change
    list(c(dR,dX,dY))
  }) #end with statement
  
}