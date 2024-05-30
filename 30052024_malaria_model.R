diff_malaria_model <- function(t, pop, para){
  
  # Assign the population matrix into the classes
  S_h <- pop[1]
  I_h <- pop[2]
  R_h <- pop[3]
  S_v <- pop[4]
  E_v <- pop[5]
  I_v <- pop[6]
  DiseaseDeaths <- 0 
  Treatment <- 0
  
  # Write down the ODE system
  dS_h <- para$mu_h*(S_h + I_h + R_h) + para$gamma*para$delta*I_h - para$a*para$p_h*S_h*I_v/(S_h + I_h + R_h) + para$omega*R_h - para$mu_h*S_h
  dI_h <- para$a*para$p_h*S_h*I_v/(S_h + I_h + R_h) - (para$gamma + para$mu_h)*I_h
  dR_h <- (1 - para$delta)*para$gamma*I_h - (para$omega + para$mu_h)*R_h
  dS_v <- para$mu_v*para$K*(1 + 0.5*cos(2*pi*t/365)) - para$a*para$p_v*S_v*I_h/(S_h + I_h + R_h) - para$mu_v*S_v
  dE_v <- para$a*para$p_v*S_v*I_h/(S_h + I_h + R_h) - (para$sigma + para$mu_v)*E_v
  dI_v <- para$sigma*E_v - para$mu_v*I_v
  
  # Return the derivatives together in a list 
  return(list(c(dS_h, dI_h, dR_h, dS_v, dE_v, dI_v)))
}

ODE_malaria_model <- function(para, ICs, maxtime){
  
  # Time points for simulation
  t_seq <- seq(0, maxtime, by = 1)
  
  # Solve the ODEs using the ode function from deSolve package
  result <- ode(y = ICs, times = t_seq, func = diff_malaria_model, parms = para, method = "ode45")
  
  # Convert the result to a data frame
  Classes <- data.frame(result)
  
  return(Classes)
}