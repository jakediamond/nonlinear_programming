# ----------------------------------------------
# Author: Jake Diamond
# Purpose: Model an RLC circuit
# Date: Februrary 23, 2017
# -----------------------------------------------

library(deSolve)
model <- function(time, states, parms){

  #States
  I1 <- states[1]
  I2 <- states[2]

  #Parameters
  r <- parms[1] #resistor (ohm)
  l <- parms[2] #inductor (h)
  c <- parms[3] #capacitor (f)

  #Parameter Calcs
  alpha <- r / (2 * l)
  w0 <- 1 / sqrt(l * c)

  #Differential Equations
  dI2 <- -2 * alpha * I2  - w0^2 * I1

  return(list(c(I2, dI2)))
}

#-------------------------------------------------
  #Define parameters and initial conditions
  parms <- c(
    r = 25/8,
    l = 0.1,
    c = 0.00000001
  )
  
  yini <- c(I1 = 2, I2 = 0)
  
  #Define time scale and simulation time
  simulation_time = 0.5 #
  dt = 0.001
  times <- seq(0, simulation_time, by = dt)
  
  #Run the model
  out = ode(y = yini, times = times, 
            func = model, parms = parms)
  
  plot(out)
  