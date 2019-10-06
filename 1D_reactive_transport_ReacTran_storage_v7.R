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
  d = 1, # average depth (m)
  w = 40, # average channel width (m)
  S = 0.001, # average channel slope (m/m)
  
  # Reactive parameters
  GPP_max = 2, # maximum GPP rate (mg O2/m2/hr)
  k_gpp = 500, # PAR at which half of GPP_max is produced (umol/m2/s)
  do_sat = 9, # saturation for DO, should be T dependent (mg/L)
  k_O2 = 1, # reaeration coefficient (1/d)
  
  # Storage parameters
  A_stor_frac = 0.2, # fraction of channel area that is storage area (-)
  alpha = 0.8, # exchange coefficient (m3/hr)
  k_er = 0.4, # storage area respiration coefficient (-)
  eta_er = 0.5, #storage area respiration efficiency (-)
  
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
model <- function(time, state, parms, gpp_choice){
  with(as.list(parms), {
    # Need to set the seed to avoid crazy computation issues, and to allow rep.
    set.seed(42)
    # Total number of boxes
    N = L / dx
    
    # Unpack states
    DO <- state[1:N]
    DO_stor <- state[(N+1):(2*N)]

    # Calculate parameters in correct units (time unit is hour)
    k = k_O2 * d / 24 # reaeration velocity (m/h)
    A = w * d # rectangular river cross-sectional area (m2)
    A_s = A * A_stor_frac # transient storage area (m2)
    
    # Calculate in-stream reaction in each box (gpp and reaeration)
    # Photosynthetically reactive radiation (umol/m2/s)
    par = ifelse(-500+2000*sin((2*pi*(time-(5/del_t)) / (24/del_t))) < 0,
                 0,
                 -500+2000*sin((2*pi*(time-(5/del_t)) / (24/del_t))))
    
    # GPP is switching function depending on user choice (g O2/m2/h)
    # First calculate base rate of GPP as Michaelis-Menten
    gpp_base = GPP_max * par / (k_gpp + par)
    # Then modulate this based on user choice
    gpp = switch(gpp_choice,
                 constant = rep(gpp_base, N),
                 ran_norm = rnorm(N, gpp_base, sd),
                 ran_uni = runif(N, 
                                 min = 0, 
                                 max = gpp_base),
                 ran_bin = gpp_base * (rbinom(N, 1, p_bi)),
                 ramp_up = seq(0, 
                               gpp_base, 
                               length.out = N),
                 ramp_down = seq(gpp_base,
                                 0, 
                                 length.out = N),
                 ramp_up_down = c(seq(0, 
                                      gpp_base, 
                                      length.out = N/2),
                                  seq(gpp_base,
                                      0, 
                                      length.out = N/2))
                 )
    # Reaeration in each box (g O2/m2/h)
    reaeration = k * (do_sat - DO)

    # advection-dispersion in each box
    adv_dis = tran.1D(C = DO,
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
    dDO = (adv_dis + (gpp + reaeration + storage_str) /d) * del_t
    dDO_stor = storage * del_t
    diffs = c(dDO = dDO, dDO_stor = dDO_stor)
    
    # Fluxes
    fluxes = c(GPP = gpp * del_t, ER = er * del_t)
    
    # Output of model
    return(list(diffs, fluxes))
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations
DO_ini = rep(7, 50)
DO_stor_ini = rep(5, 50)
yini <- c(DO = DO_ini, DO_stor = DO_stor_ini)

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 0.25 # time step (h)
days <- 4 # number of days to simulate
simulation_time <- 24 * days / del_t # simulation time (h)
times <- seq(0, simulation_time, by = del_t)

# GPP choice (choose between: "constant", "ran_uni", "ran_norm", "ran_bin", 
# "ramp_up", "ramp_down", or "ramp_up_down")
gpp_choice <- "ran_uni"

# uses the R deSolve function (lsoda method)
out <- ode.1D(y = yini, 
              times = times,
              func = model, 
              parms = parms,
              nspec = 2,
              dimens = with(as.list(parms), L / dx),
              gpp_choice = gpp_choice)

# Examine the model output ------------------------------------------------
# Reorganize the data into long-form
df <- as_tibble(out) %>%
  gather(key, value, -time) %>%
  separate(key, c("type", "reach"), "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(time_hr = time * del_t,
         dist = as.numeric(reach) * with(as.list(parms), dx)) %>%
  mutate(type_plot = recode(type,
         `DO` = "DO~(mg~L^{-1})",
         `DO_stor` = "DO[storage]~(mg~L^{-1})",
         `ER.DO_stor` = "ER~(g~O^2~m^{-2}~h^{-1})",
         `GPP` = "GPP~(g~O^2~m^{-2}~h^{-1})"))

# Plot the data
ggplot(data = filter(df, reach %in% c(10, 20, 30, 40)),
       aes(x = time_hr,
           y = value,
           color = dist)) +
  geom_point(alpha = 0.4) +
  facet_grid(rows = vars(type_plot), scales = "free_y",
             labeller = label_parsed) + 
  scale_color_viridis_c(name = "Distance (m)") +
  scale_x_continuous(breaks = seq(0, simulation_time * del_t, 24)) + 
  theme_bw() +
  xlab("Time (h)") +
  ylab("")

# Summarize the data
df %>%
  group_by(type) %>%
  summarize(avg = mean(value))

