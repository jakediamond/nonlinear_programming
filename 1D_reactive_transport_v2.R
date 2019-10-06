# 1D advection, dispersion, reaction model for dissolved oxygen 

library(deSolve)
library(tidyverse)

model <- function(time, state, parms){
  # Unpack states
  DO = state[1:N]
  
  # Calculate reaction in each box (gpp, er, and reaeration)
  par = ifelse(540+800*sin(2*pi*time/24-1.4) < 0, # model PAR ()
               0,
               540+800*sin(2*pi*time/24-1.4))
  reaeration = k_O2 * (do_sat - DO) # model reaeration
  gpp = k_gpp * par #/ d / 1000 # model GPP (g O2/m2/d)
  er = -0.1#-0.5 * gpp # model ER (g O2/m2/d)
  reaction = reaeration + gpp + er
  
  # Advection
  # deltax <- c (0.5, rep(1, N-2), 0.5)
  advection = if(time == 0){
    rep(0, N)
  } else {
    -(Q / A) * 3600 * (lead(DO, default = 0) - 
                         lag(DO, default = 0)) / (2 * dx)
  }
  
  # Dispersion
  dispersion = if(time == 0){
    rep(0, N)
  } else {
    (D * (lead(DO) - (2 * DO) + lag(DO))) / dx^2
  }
  
  # Overall difference
  dDO = advection + dispersion + reaction
  
  return(list(c(dDO = dDO)))     # result
}  # end of model


# Model parameters --------------------------------------------------------

# Grid and time parameters
dx <- 80 # segment length (m)
L <- 4000 # reach length (m)
# x <- seq(dx / 2, L, by = dx) # distance along reach (m)
N <- L / dx # number of segments

# Hydraulic parameters
Q <- 5 # discharge (m3/s)
d <- 1.5 # average depth (m)
w <- 40 # average channel width (m)
S <- 0.001 # average channel slope (m/m)
A = w * d # rectangular river cross-sectional area (m2)
u_s = sqrt(9.8 * d * S) # shear velocity
D = 4#0.6 * d *  u_s # longitudinal dispersion (m2/s)

# Reactive parameters
k_gpp <- 0.00051 # coefficient for converting PAR to GPP (mg O2/J)
do_sat <- 9 # saturation for DO, should be T dependent (mg/L)
k_O2 <- 1 # reaeration coefficient (1/d)

# Initial conditions ------------------------------------------------------
state <- rep(do_sat - 1, N) 
  DO = rnorm(N, do_sat, 1) # set initial DO to just below saturation (mg/L)

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 1/6 # time step (hr)
simulation_time <- 24 * 2 # hours
times <- seq(0, simulation_time, by = del_t)

# uses the R deSolve function (adaptive Runge-Kutta method)
out <- ode.1D(y = state, 
           times = times, 
           func = model, 
           parms = NULL, 
           nspec = 1,
           names = "DO",
           dimens = N)

# Examine the model output ------------------------------------------------
# Time series of every 5 segments
as_tibble(out) %>%
  select(time, seq(2, ncol(out), 4)) %>%
  gather(location, value, -time) %>%
  ggplot(aes(x = time, y = value, color = location)) +
  geom_point()

