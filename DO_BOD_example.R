O2BOD <- function(t, state, pars) {
  BOD <- state[1:N]
  O2 <- state[(N+1):(2*N)]
  ## BOD dynamics
  FluxBOD <- v * c(BOD_0, BOD) # fluxes due to water transport
  FluxO2 <- v * c(O2_0, O2)
  BODrate <- r * BOD # 1-st order consumption
  ## rate of change = flux gradient - consumption + reaeration (O2)
  dBOD <- -diff(FluxBOD)/dx - BODrate
  dO2 <- -diff(FluxO2)/dx - BODrate + p * (O2sat-O2)
  return(list(c(dBOD = dBOD, dO2 = dO2)))
}


dx <- 25 # grid size of 25 meters
v <- 1e3 # velocity, m/day
x <- seq(dx/2, 5000, by = dx) # m, distance from river
N <- length(x)
r <- 0.05 # /day, first-order decay of BOD
p <- 0.5 # /day, air-sea exchange rate
O2sat <- 300 # mmol/m3 saturated oxygen conc
O2_0 <- 200 # mmol/m3 riverine oxygen conc
BOD_0 <- 1000 # mmol/m3 riverine BOD concentration

## initial conditions:
state <- c(rep(200, N), rep(200, N))
times <- seq(0, 20, by = 0.1)
## running the model
## step 1 : model spinup
out <- ode.1D(y = state, times, O2BOD, parms = NULL,
              nspec = 2, names = c("BOD", "O2"))

## select oxygen (first column of out:time, then BOD, then O2
O2 <- out[, (N + 2):(2 * N + 1)]
color = topo.colors
filled.contour(x = times, y = x, O2, color = color, nlevels = 50,
               xlab = "time, days", ylab = "Distance from river, m",
               main = "Oxygen")
## or quicker plotting:
image(out, grid = x, xlab = "time, days", ylab = "Distance from river, m")
