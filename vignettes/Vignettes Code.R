# Vignettes Code
# instructions on obtaining `fftw`



# Example of Simulation and Estimation ------------------------------------
setwd("D:/GitHub/SuperGauss/vignettes")
require(numDeriv)
require(SuperGauss)
require(fftw)
source("acf-functions.R")

# Simulating Data
# model: fBM + stochastic drift
N <- 2000
d <- 5
dT <- 1/60
alpha <- 0.4
lambda <- 100 * dT
mu <- 0.01
acf.fBM <- fbm.acf(alpha, dT, N)
acf.exp <- exp2.acf(lambda, dT, N)
fBM <- t(rSnorm(d, acf = acf.fBM))
exp <- t(rSnorm(d, acf = acf.exp))
X <- mu * dT * 1:N + fBM + exp

# Estimating Parameters
# using Newton-Raphson



# Simulating Posterior Distribution
# using Dietrich algorithm