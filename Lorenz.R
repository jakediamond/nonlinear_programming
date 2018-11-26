


dx <- sigma * (y - x)
dy <- rho * x - y - x * z
dz <- x * y - beta * z

sigma <- 10
rho <- 28
beta <- 8/3



#######matrix form

dy <- A * y
y <- c(y1, y2, y3)
A <- matrix(data = c(-beta, 0, y2, 
                     0, -sigma, sigma,
                     -y2, rho, -1), nro = 3, ncol =3)

# Make A singular
eta <- sqrt(beta * (rho - 1))
eta <- -sqrt(beta * (rho - 1))
Y.0 <- c(rho - 1, eta, eta)

A <- matrix(data = c(-beta, 0, eta, 
                     0, -sigma, sigma,
                     -eta, rho, -1), nro = 3, ncol =3)



parameters <- c(s = 10, r = 600, b = 8/3)
state <- c(X = 0, Y = 1, Z = 1)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y
    dZ <- X * Y - b * Z
    list(c(dX, dY, dZ))
  })
}
times <- seq(0, 50, by = 0.01)
library(deSolve)
out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

par(oma = c(0, 0, 3, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "Y"], out[, "Z"], pch = ".", type = "l")
mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)
