#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To model dissolved oxygen spatiotemporal evolution in streams/rivers
# Date: October 5, 2019
# 

# Load libraries
library(deSolve)
library(tidyverse)
library(ReacTran)

# Define model parameters --------------------------------------------------------
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
  k_gpp = 0.00051, # coefficient for converting PAR to GPP (mg O2/L/J)
  do_sat = 9, # saturation for DO, should be T dependent (mg/L)
  k_O2 = 1, # reaeration coefficient (1/d)
  
  # Storage parameters
  A_stor_frac = 0.2, # fraction of channel area that is storage area (-)
  alpha = 0.278, # exchange coefficient (m3/hr)
  k_er = 0.000398, # storage area respiration coefficient (-)
  eta_er = 0.02, #storage area respiration efficiency (-)
  
  # Random GPP parameters
  sd = 0.2, # standard deviation of normal distribution to GPP
  p_bi = 0.5 # probability for binomial distribution of GPP
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
model <- function(time, yini, parms, gpp_choice){
  with(as.list(parms), {
    # Need to set the seed to avoid crazy computation issues, and to allow rep.
    set.seed(42)
    # Total number of boxes
    N = L / dx
    # Unpack states
    DO = yini[1:N]
    DO_stor = yini[(N+1):(2*N)]
    
    # Calculate parameters in correct units (time unit is hour)
    k = k_O2 / 24 # reaeration velocity (1/h)
    A = w * d # rectangular river cross-sectional area (m2)
    A_s = A * A_stor_frac # transient storage area (m2)
    
    # Calculate in-stream reaction in each box (gpp and reaeration)
    # Photosynthetically reactive radiation (J/h)
    par = ifelse(540+800*sin(2*pi*time/24-1.4) < 0,
                 0,
                 540+800*sin(2*pi*time/24-1.4))
    # GPP is switching function depending on user choice (mg O2/L/h)
    gpp = switch(gpp_choice,
                 constant = k_gpp * par,
                 ran_norm = rnorm(N, k_gpp * par, sd),
                 ran_uni = runif(N, min = 0, max = k_gpp * par),
                 ran_bin = k_gpp * par * (rbinom(N, 1, p_bi)),
                 ramp_up = seq(0, k_gpp * par, by = (k_gpp * par / (N-1))),
                 ramp_down = seq(k_gpp * par, 0, by = -(k_gpp * par / (N-1))),
                 parab_up = c(seq(0, 
                               sqrt(k_gpp * par), 
                               by = (sqrt(k_gpp * par)) / (N / 2 - 1)) ^ 2,
                           seq(sqrt(k_gpp * par), 
                               0, 
                               by = -(sqrt(k_gpp * par)) / (N / 2 -1)) ^ 2
                           ),
                 parab_down = c(seq(sqrt(k_gpp * par), 
                                    0, 
                                    by = -(sqrt(k_gpp * par)) / (N / 2 -1)) ^ 2,
                                seq(0, 
                                    sqrt(k_gpp * par), 
                                    by = (sqrt(k_gpp * par)) / (N / 2 - 1)) ^ 2
                                )
                 )
    reaeration = k * (do_sat - DO) # (mg O2/L/h)
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
    storage = alpha * (A / A_s) * (DO - DO_stor) - er
      
    # Change in stream DO concentration due to transient storage
    storage_str = alpha * (DO_stor - DO)
    
    # Rates of change
    dDO = advection + reaction + storage_str
    dDO_stor = storage

    return(list(c(dDO = dDO, dDO_stor = dDO_stor),
                c(GPP = gpp, ER = er)))  # results
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations
yini <- c(DO = seq(7, 9, (2 / (50 - 1))),
             DO_stor = seq(7, 9, (2 / (50 - 1))))

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 10/60 # time step (h)
days <- 7 # number of days to simulate
simulation_time <- 24 * days # simulation time (h)
times <- seq(0, simulation_time, by = del_t)

# GPP choice (choose between: "constant", "ran_uni", "ran_norm", "ramp_up", 
# "ramp_down", or "parab")
gpp_choice = "parab_down"

# uses the R deSolve function (lsoda method)
out <- ode.1D(y = yini, 
              times = times,
              func = model, 
              parms = parms,
              nspec = 2,
              gpp_choice = gpp_choice)

# Examine the model output ------------------------------------------------
# Time series of every 5 segments 
# If you choose "constant" as the gpp_choice, they are
# the data in 'out' consist of: 1 = times; 2:N+1 = DO; N+2:2*N+1 = DO_stor; 
# 2*N+2 = GPP; and 2*N+3:3*N+1 = ER
as_tibble(out) %>%
  # select(time, seq(2, ncol(out) / 2, 20)) %>% # DO in stream at reach 1, 20, 40
  select(time, 102) %>% #GPP
  gather(location, value, -time) %>%
  ggplot(aes(x = time, y = value, color = location)) +
  geom_point()

