# 1D advection, dispersion, reaction model for dissolved oxygen 

library(deSolve)
library(tidyverse)
library(ReacTran)

# Model parameters --------------------------------------------------------
parms <- c(
  # Grid and time parameters
  dx = 80, # segment length (m)
  L = 4000, # reach length (m)
  
  # Hydraulic parameters
  Q = 5, # discharge (m3/s)
  d = 1.5, # average depth (m)
  w = 40, # average channel width (m)
  S = 0.001, # average channel slope (m/m)
  
  # Reactive parameters
  k_gpp = 0.00051, # coefficient for converting PAR to GPP (mg O2/J)
  do_sat = 9, # saturation for DO, should be T dependent (mg/L)
  k_O2 = 1, # reaeration coefficient (1/d)
  
  # Storage parameters
  A_stor_frac = 0.2, # fraction of channel area that is storage area (-)
  alpha = 0.278, # exchange coefficient (m3/hr)
  k_er = 0.000398, # storage area respiration coefficient (-)
  eta_er = 0.02 #storage area respiration efficiency (-)
 )

# Set up model grid -------------------------------------------------------
# Grid of reaches with total length L and reach length L / dx
grid <- with(as.list(parms),
             setup.grid.1D(L = L, N = L / dx)
             )
# Dispersion grid, all the same longitudinal dispersion based on shear velocity
D.grid   <- with(as.list(parms),
                 setup.prop.1D(value = d * sqrt(9.8 * d * S) * 3600, # (m2/h)
                               grid = grid)
                 ) 
# Velocity grid, all the same based on discharge and area (assumed rectangular)
v.grid   <- with(as.list(parms),
                 setup.prop.1D(value = Q * 3600 / (w * d), # (m/hr)
                               grid = grid)
                 )

# Overall model function ----------------------------------------------------------
model <- function(time, yini, parms){
  with(as.list(parms), {
    # Total number of boxes
    N = L / dx
    # Unpack states
    DO = yini[1:N]
    DO_stor = yini[(N+1):(2*N)]
    
    # Get parameters in correct units
    k = k_O2 * d / 24 # reaeration velocity (m/h)
    A = w * d # rectangular river cross-sectional area (m2)
    A_s = A * A_stor_frac # transient storage area (m2)
    
    # Calculate reaction in each box (gpp, er, and reaeration)
    par = ifelse(540+800*sin(2*pi*time/24-1.4) < 0,
                 0,
                 540+800*sin(2*pi*time/24-1.4))
    reaeration = k * (do_sat - DO)
    gpp = k_gpp * par # g O2/m2/d
    reaction = reaeration + gpp
    
    # advection-dispersion in each box
    advection = tran.1D(C = DO,
                        C.up = do_sat, # upstream boundary concentration
                        C.down = do_sat - 1, # downstream boundary concentration
                        D = D.grid, 
                        v = v.grid, 
                        dx = grid)$dC 
    
    # Change in storage concentration due to respiration,
    # which only happens in transient storage 
    er = k_er * DO_stor ^ eta_er
    storage = alpha * A * (DO - DO_stor) / A_s - er
      
    # Change in stream DO concentration due to transient storage
    storage_str = alpha * (DO_stor - DO)
    
    # Rates of change
    dDO = advection + reaction + storage_str
    dDO_stor = storage

    return(list(c(dDO = dDO, dDO = dDO_stor),
                c(GPP = gpp, ER = er)))  # results
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations
yini <- c(DO = seq(7, 9, (2 / (50 - 1))),
             DO_stor = seq(7, 9, (2 / (50 - 1))))

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 10/60 # time step (hr)
simulation_time <- 24 * 7 # hours
times <- seq(0, simulation_time, by = del_t)

# uses the R deSolve function (adaptive Runge-Kutta method)
out <- ode.1D(y = yini, 
              times = times,
              func = model, 
              parms = parms,
              nspec = 2)

# Examine the model output ------------------------------------------------
# Time series of every 5 segments
# # the data in 'out' consist of: 1st col times, 2:N+1 DO;
# N+2:2*N+1 GPP; and 2*N+2:3*N+1 ER
as_tibble(out) %>%
  # rename(2:51,) %>%
  select(time, seq(2, ncol(out) / 2, 20)) %>%
  gather(location, value, -time) %>%
  ggplot(aes(x = time, y = value, color = location)) +
  geom_point()

