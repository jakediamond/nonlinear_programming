# 1D advection, dispersion, reaction model for dissolved oxygen 

library(deSolve)
library(tidyverse)
library(ReacTran)

# Model parameters --------------------------------------------------------
parms <- c(
  # Grid and time parameters
  dx <- 80, # segment length (m)
  L <- 4000, # reach length (m)
  # x <- seq(dx / 2, L, by = dx) # distance along reach (m)
  N <- L / dx, # number of segments
  
  # Hydraulic parameters
  Q <- 5, # discharge (m3/s)
  d <- 1.5, # average depth (m)
  w <- 40, # average channel width (m)
  S <- 0.001, # average channel slope (m/m)
  A <- w * d, # rectangular river cross-sectional area (m2)
  u_s <- sqrt(9.8 * d * S), # shear velocity
  D <- d * u_s, # longitudinal dispersion (m2/s)
  
  # Reactive parameters
  k_gpp <- 0.00051, # coefficient for converting PAR to GPP (mg O2/J)
  do_sat <- 9, # saturation for DO, should be T dependent (mg/L)
  k_O2 <- 0.00000694471905785 # reaeration coefficient (1/d)
)


# Set up model grid -------------------------------------------------------
grid <- setup.grid.1D(L = L, N = N)
D.grid   <- setup.prop.1D(value = D, grid = grid)
v.grid   <- setup.prop.1D(value = Q/A, grid = grid)

# Model function ----------------------------------------------------------
model <- function(time, DO, parms){
  with (as.list(parms), {
    # Calculate reaction in each box (gpp, er, and reaeration)
    par = ifelse(540+800*sin(2*pi*time/24-1.4) < 0,
                 0,
                 540+800*sin(2*pi*time/24-1.4))
    reaeration = k_O2 * (do_sat - DO)
    gpp = k_gpp * par + rnorm(length(DO), 0, 0.3) #/ d / 1000 #g O2/m2/d
    er = -0.1 + -0.5 * gpp
    reaction = reaeration + gpp + er
    
    # advection-dispersion in each box
    advection = tran.1D(C = DO,
                        C.up = do_sat,
                        C.down = do_sat + 1,
                        D = D.grid, 
                        # F.up = NULL, 
                        # A = 1,
                        v = v.grid, 
                        dx = grid)$dC * 1000
    
    # Rate of change
    dDO = advection + reaction

    return(list(dDO = dDO, GPP = gpp, ER = er))     # result
  })
}  # end of model

# Initial conditions ------------------------------------------------------
yini <- seq(do_sat - 2, do_sat, (2 / (N - 1)))

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 0.5 # time step (hr)
simulation_time <- 24 * 7 # hours
times <- seq(0, simulation_time, by = del_t)

# uses the R deSolve function (adaptive Runge-Kutta method)
out <- ode.1D(y = yini, 
           times = times,
           func = model, 
           parms = parms,
           nspec = 1,
           method = "ode45")

# Examine the model output ------------------------------------------------
# Time series of every 5 segments
# # the data in 'out' consist of: 1st col times, 2:N+1 DO;
# N+2:2*N+1 GPP; and 2*N+2:3*N+1 ER
as_tibble(out) %>%
  select(time, seq(2, ncol(out) / 3, 20)) %>%
  gather(location, value, -time) %>%
  ggplot(aes(x = time, y = value, color = location)) +
  geom_point()



