#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To model dissolved oxygen spatiotemporal evolution in streams/rivers
# Date: October 5, 2019
# 

# Load libraries
library(deSolve)
library(tidyverse)
library(ReacTran)
library(lubridate)

# Define model parameters --------------------------------------------------------
parms <- c(
  # Grid and time parameters
  dx = 80, # segment length (m)
  L = 4000, # reach length (m)
  C_down = 9, # downstream boundary concentration (mg/L)
  C_up = 9, # upstream boundary concentration (mg/L)
  
  # Hydraulic parameters
  Q = 5, # discharge (m3/s)
  d = 1.5, # average depth (m)
  w = 40, # average channel width (m)
  S = 0.001, # average channel slope (m/m)
  
  # Reactive parameters
  GPP_max = 1.5, # maximum GPP rate (g O2/m2/h)
  k_gpp = 500, # PAR at which half of GPP_max is produced (umol/m2/s)
  do_sat = 9, # saturation for DO, should be T dependent (mg/L)
  k_O2 = 0.6, # reaeration coefficient (1/d)
  er_const = -0.393, # constant ER rate (g O2/m2/h)
  chan_frac = 0.5, # fraction of ER that occurs in channel vs storage (-)
  
  # Storage parameters
  A_stor_frac = 0.2, # fraction of channel area that is storage area (-)
  alpha = 0.1, # exchange coefficient (1/h)
  k_er = 0.4, # storage area respiration coefficient (-)
  eta_er = 0.2, #storage area respiration efficiency (-)
  
  # Random GPP parameters
  sd = 0.2, # standard deviation of normal distribution to GPP
  p_bi = 0.5 # probability for binomial distribution of GPP
 )

# Set up model grid -------------------------------------------------------
# Grid of reaches with total length L and reach length L / dx
grid <- with(as.list(parms),
             setup.grid.1D(L = L, N = L / dx)
             )
# Dispersion grid, all the same longitudinal dispersion based on shear velocity (Kashefipour and Falconer 2002)
D.grid   <- with(as.list(parms),
                 setup.prop.1D(value = 14400, # (m2/h)
                               grid = grid)
                 ) 

# Velocity grid, all the same based on discharge and area (assumed rectangular)
v.grid   <- with(as.list(parms),
                 setup.prop.1D(value = Q * 3600 / (w * d), # (m/h)
                               grid = grid)
                 )
# 10.612 * (Q / (w * d))^2 / (sqrt(9.8 * d * S)) * 3600
# Overall model function ----------------------------------------------------------
model <- function(time, state, parms, gpp_choice, er_loc_choice, er_choice){
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
    # Photosynthetically active radiation (umol/m2/s)
    par = ifelse(-500+2000*sin((2*pi*(time-(6/del_t)) / (24/del_t))) < 0,
                 0,
                 -500+2000*sin((2*pi*(time-(6/del_t)) / (24/del_t))))
    
    # GPP is switching function depending on user choice (g O2/m2/h)
    # First calculate base rate of GPP as Michaelis-Menten
    gpp_base = GPP_max * par / (k_gpp + par)
    # Then modulate this based on user choice
    gpp = switch(gpp_choice,
                 none = rep(0, N),
                 constant = rep(gpp_base, N),
                 ran_norm = rnorm(N, gpp_base, sd),
                 ran_uni = runif(N, 
                                 min = 0.5*gpp_base, 
                                 max = 1.5*gpp_base),
                 ran_bin = gpp_base * (rbinom(N, 1, p_bi)),
                 ramp_up = seq(0, 
                               gpp_base*2, 
                               length.out = N),
                 ramp_down = seq(gpp_base*2,
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

    # advection-dispersion in each box (g O2/m3/h)
    adv_dis = tran.1D(C = DO,
                      C.up = C_up, # upstream boundary concentration
                      C.down = C_down, # downstream boundary concentration
                      D = D.grid, 
                      v = v.grid, 
                      dx = grid)$dC 
    
    # Ecosystem respiration amount (g O2/m2/h)
    er = switch(er_choice,
                none = rep(0, N),
                constant = rep(er_const, N),
                efficiency = -k_er * DO_stor ^ eta_er
                )
    
    # Ecosystem respiration location changes where DO is consumed
    # Change in storage concentration due to respiration (g O2/m3/h)
    storage = switch(er_loc_choice,
                     storage = alpha * (A / A_s) * (DO - DO_stor) + 
                       er / (A_s / w),
                     channel = alpha * (A / A_s) * (DO - DO_stor),
                     both = alpha * (A / A_s) * (DO - DO_stor) + 
                       (er * (1 - chan_frac) / (A_s / w))
                     )
    
    # Change in stream DO concentration due to transfer to transient storage (g O2/m3/h)
    storage_str = alpha * (DO_stor - DO)
    
    # Rates of change
    dDO = switch(er_loc_choice,
                 storage = (adv_dis + 
                              storage_str + (gpp + reaeration ) / d) * del_t,
                 channel = (adv_dis + 
                              storage_str + (gpp + reaeration + er) / d) * del_t,
                 both = (adv_dis + storage_str + 
                           (gpp + reaeration + 
                              er * chan_frac
                            ) / d
                         ) * del_t,
                 )
    dDO_stor = storage * del_t
    diffs = c(dDO = dDO, dDO_stor = dDO_stor)
    
    # Fluxes
    fluxes = c(GPP = gpp * del_t, ER = er * del_t)
    
    # Output of model
    return(list(diffs, fluxes, PAR = rep(par, N)))
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations
DO_ini = 9 # mg/L
DO_stor_ini = 9 # mg/L
yini <- c(DO = rep(DO_ini, with(as.list(parms), L / dx)),
          DO_stor = rep(DO_stor_ini, with(as.list(parms), L / dx)))

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 1 # time step (h)
days <- 8 # number of days to simulate
simulation_time <- 24 * days / del_t # simulation time (h)
times <- seq(0, simulation_time, by = del_t)

# GPP choice (choose between: "none", "constant", "ran_uni", "ran_norm", "ran_bin", 
# "ramp_up", "ramp_down", or "ramp_up_down")
gpp_choice <- "constant"

# ER location choice (choose between: "storage", "channel", or "both")
er_loc_choice <- "storage"

# ER choice (choose between: "constant", or "efficiency")
er_choice <- "constant"

# uses the R deSolve function (lsoda method)
out <- ode.1D(y = yini, 
              times = times,
              func = model, 
              parms = parms,
              nspec = 2,
              dimens = with(as.list(parms), L / dx),
              gpp_choice = gpp_choice,
              er_loc_choice = er_loc_choice,
              er_choice = er_choice)

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
         `ER` = "ER~(g~O^2~m^{-2}~h^{-1})",
         `GPP` = "GPP~(g~O^2~m^{-2}~h^{-1})",
         `PAR` = "PAR~({`mu`}*mol~m^{-2}~s^{-1})"))

# Plot the data
ggplot(data = filter(df, reach %in% c(20, 40, 60)),
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
  summarize(avg = mean(value) * 24)

# Get data in correct form for streammetabolizer --------------------------
df_sm <- df %>%
  filter(type %in% c("PAR", "DO"),
         reach == 95) %>%
  select(time_hr, type, value) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  rename(light = PAR,
         DO.obs = DO) %>%
  mutate(temp.water = 20,
         DO.sat = with(as.list(parms), do_sat),
         depth = with(as.list(parms), d),
         day_no = floor(time_hr / 24) + 1,
         time = paste(floor(time_hr %% 24), 
                      round((time_hr %% 24 - floor(time_hr %% 24)) * 60,
                            digits = 2), sep=":"),
         datetime = ymd_hm(paste0("2019-07-0", day_no, " ", time))) %>%
  select(-time_hr, -time, -day_no)




x <- df %>%
  filter(type %in% c("GPP", "ER"),
         reach == 95) %>%
  select(time_hr, type, value) %>%
  spread(type, value) %>%
  mutate(
         day_no = floor(time_hr / 24) + 1,
         time = paste(floor(time_hr %% 24), 
                      round((time_hr %% 24 - floor(time_hr %% 24)) * 60,
                            digits = 2), sep=":"),
         datetime = ymd_hm(paste0("2019-07-0", day_no, " ", time)))
x <- x %>% 
  group_by(day_no) %>%
  summarize(gpp = sum(GPP),
            er = sum(ER))
