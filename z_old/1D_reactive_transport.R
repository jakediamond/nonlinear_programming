# 1D advection, dispersion, reaction model for dissolved oxygen 

library(deSolve)

model <- function(time, states, params){
  # Unpack states
  DO = states
  dDO = rep(NA,length(DO))
  # Unpack parameters
  N = params[1]
  L = params[2]
  Q = params[3]
  d = params[4]
  w = params[5]
  S = params[6]
  k_gpp = params[7]
  do_sat = params[8]
  k_O2 = params[9]
  
  # Calculate other parameters
  del_x = L/N # segment length (m)
  A = w * d # rectangular river cross-sectional area (m2)
  u_s = sqrt(9.8 * d * S) # shear velocity
  D = 0.6 * d *  u_s # longitudinal dispersion (m2/s)
  
  # Loop through segments
  for(i in 1:(N)){
    
    # Calculate reaction in each box (gpp, er, and reaeration)
    par = ifelse(540+800*sin(2*pi*time/24-1.4) < 0, # model PAR ()
                 0,
                 540+800*sin(2*pi*time/24-1.4))
    reaeration = k_O2 * (do_sat - DO[i]) # model reaeration
    gpp = k_gpp * par #/ d / 1000 # model GPP (g O2/m2/d)
    er = -0.1#-0.5 * gpp # model ER (g O2/m2/d)
    reaction = reaeration + gpp + er

    if(i == 1){  # upper boundary condition
      # Calculate upper boundary layer advection
      advection = 0
      # Calculate upper boundary layer dispersion
      dispersion = 0
    } else if(i > 1 & i < N){  # boxes in the middle
      # Calculate mid reach advection
      advection = -(Q / A) * 3600 * (lead(DO, default = 0) - 
                            lag(DO, default = 0)) / (2 * dx)
      # Calculate mid reach dispersion
      dispersion = (D * (lead(DO) - (2 * DO) + lag(DO))) / dx^2
    } else{  # lower boundary condition
      # Calculate lower boundary advection
      advection = -(Q / A) * dC / del_x
      # Calculate lower boundary dispersion
      dispersion = eppd
    }
    
    # Calculate derivatives
    dDO[i] = reaction + advection + dispersion
  }
  
  list(dDO)     # result
  
}  # end of model


# Model parameters --------------------------------------------------------
params <- c(
  # Grid and time parameters
  N <- 10, # number of segments
  L <- 100, # reach length (m)

  # Hydraulic parameters
  Q <- 0.1, # discharge (m3/s)
  d <- 0.5, # average depth (m)
  w <- 10, # average channel width (m)
  S <- 0.001, # average channel slope (m/m)

  # Reactive parameters
  k_gpp <- 0.00051, # coefficient for converting PAR to GPP (mg O2/J)
  do_sat <- 9, # saturation for DO, should be T dependent (mg/L)
  k_O2 <- 0.1 # reaeration coefficient (1/d)
)
# Initial conditions ------------------------------------------------------
yini <- c(
  DO = rep(8, times = N) # set initial DO to just below saturation (mg/L)
  )
# # state <- rep(0, times = N)
# state       <- rep(0,times=numboxes)
# 
# state[1]=1
# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 0.5 # time step (hr)
simulation_time <- 24 * 7 # hours
times <- seq(1, simulation_time, by = del_t)

# uses the R deSolve function (adaptive Runge-Kutta method)
out <- ode(y = yini, 
           times = times, 
           func = model, 
           parms = params, 
           method = "ode45")

# Examine the model output ------------------------------------------------
# Time series of every 5 segments
as_tibble(out) %>%
  select(time, seq(2, ncol(out), 4)) %>%
  gather(location, value, -time) %>%
  ggplot(aes(x = time, y = value, color = location)) +
  geom_point()



# the data in 'out' consist of: 1st col times, 2-41: the concentrations
# Select the concentration data
result <- out[, 2:(N + 1)]

length  <- seq(from = 0, by = del_x, length.out = N)
# 1. temporal-spatial plot of the concentrations
par(oma=c(0,0,3,0))   # set margin size (oma)
col <- topo.colors
#col <-   greycol
filled.contour(x = times,
               y = length,
               result,
               color = col,
               # ylim = c(9.5, 0.5),
               # zlim = c(6, 12),
               xlab = "time, days", 
               ylab = "Depth, m",
               main = "Concentration, gC/m3")

mtext(outer=TRUE,side=3,"Vertical Phyto model",cex=1.5)

phyto_total = rowSums(PHYTO)

plot(times,phyto_total,type='l',xlab='times, days',ylab='Concentration, gC/m3')

