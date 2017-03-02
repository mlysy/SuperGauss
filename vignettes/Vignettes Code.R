# Vignettes Code

# Instruction On FFTW -----------------------------------------------------



# Example of Simulation and Estimation ------------------------------------

setwd("D:/GitHub/SuperGauss/vignettes")
require(numDeriv)
require(SuperGauss)
require(fftw)
source("acf-functions.R")

# Data Simulation ---------------------------------------------------------

# model: fBM + stochastic drift
N <- 2000
d <- 1
dT <- 1/60
H <- 0.4
lambda <- 200 * dT
mu <- 0.5
theta <- c(mu, H, lambda)

d.drift <- mean.fun(mu, dT, N)
d.fbm <- rSnorm(d, acf = fbm.acf(H, dT, N))
d.exp2 <- rSnorm(d, acf = exp2.acf(lambda, dT, N))

dX <- matrix(d.drift + d.fbm + d.exp2, N)
drift <- cumsum(matrix(d.drift + d.fbm, N))
vol <- cumsum(matrix(d.exp2, N))
X <- drift + vol
rng <- range(c(X, drift, vol))
plot(X, col = "blue", type = "l", ylim = rng, xlab = "", ylab = "")
par(new = TRUE)
plot(drift, col = "red", type = "l", ylim = rng, xlab = "", ylab = "")
par(new = TRUE)
plot(vol, col = "green", type = "l", ylim = rng, xlab = "", ylab = "")

# Parameter Estimation Using Newton-Raphson -------------------------------

acf.fun <- function(H, lambda, dT = dT, N = N){
  acf1 <- fbm.acf(H, dT, N)
  acf2 <- exp2.acf(lambda, dT, N)
  matrix(acf1 + acf2, N, 1)
}

d.mu <- function(mu, dT, N){
  2 * mu * dT * matrix(1, N, 1)
}

d.fbm <- function(H, dT = dT, N = N){
  acf <- 2 * (dT*(0:N+1))^(2*H) * log(dT*(0:N+1))
  acf <- -1/2 * acf2incr(acf)
  matrix(acf, N, 1)
}

d.exp2 <- function(lambda, dT = dT, N = N){
  acf <- 2 * (dT*(0:N))^2 * exp(-(dT*(0:N)/lambda)^2) / lambda^3
  acf <- acf2incr(acf)
  matrix(acf, N, 1)
}

d2.mu <- function(mu, dT, N){
  2 * dT * matrix(1, N, 1)
}

d2.fbm <- function(H, dT = dT, N = N){
  acf <- 4 * (dT*(0:N+1))^(2*H) * log(dT*(0:N+1))^2
  acf <- -1/2 * acf2incr(acf)
  matrix(acf, N, 1)
}

d2.exp2 <- function(lambda, dT = dT, N = N){
  acf <- (4 * (dT*(0:N))^4 - 6 * (dT*(0:N))^2 * lambda^2) * exp(-(dT*(0:N)/lambda)^2) / lambda^6
  acf <- acf2incr(acf)
  matrix(acf, N, 1)
}

# Newton-Raphson Using Function 'nlm'

dSnorm.para <- function(theta, X, dT, Toep){
  n <- nrow(X)
  mu <- theta[1]
  H <- theta[2]
  lambda <- theta[3]
  p <- length(theta)
  mean <- mean.fun(mu, dT, n)
  acf <- acf.fun(H, lambda, dT, n)
  
  dmean <- matrix(0, n, p)
  dmean[, 1] <- d.mu(mu, dT, n)
  
  dacf <- matrix(0, n, p)
  dacf[, 2] <- d.fbm(H, dT, n)
  dacf[, 3] <- d.exp2(lambda, dT, n)
  
  d2mean <- array(0, c(n, 3, 3))
  d2mean[, 1, 1] <- d2.mu(mu, dT, n)
  
  d2acf <- array(0, c(n, 3, 3))
  d2acf[, 2, 2] <- d2.fbm(H, dT, n)
  d2acf[, 3, 3] <- d2.exp2(lambda, dT, n)
  
  X1 <- X - mean
  Toep$AcfInput(acf)
  density <- crossprod(X1, Toep$Solve(X1))
  density <- (density + Toep$Det()) / 2
  attr(density, "gradient") <- -1 * Snorm.grad(X, mean, acf, dmean, dacf, Toep)
  attr(density, "hessian") <- -1 * Snorm.Hess(X, mean, acf, dmean, dacf, d2mean, d2acf, Toep)
  density
}

Toep <- new(Toeplitz, N)
theta.start <- c(0.1, 0.5, 3)
theta.nlm <- nlm(f = dSnorm.para, p = theta.start, X = dX, dT = dT, Toep = Toep)
signif(cbind(theta.nlm$estimate, theta))


# Drift Simulation --------------------------------------------------------

# Drift simulation technique with true parameter
dT <- 1/60
N <- 2000
param <- theta
mu.para <- param[1]
H.para <- param[2]
lambda.para <- param[3]
Toep <- new(Toeplitz, N)
mean <- mean.fun(mu.para, dT, N)
acf1 <- fbm.acf(H.para, dT, N)
acf2 <- exp2.acf(lambda.para, dT, N)

# simu
nSim <- 100
path.sim <- matrix(NA, N, nSim)
e1 <- rSnorm(nSim, acf = acf1)
e2 <- rSnorm(nSim, acf = acf2)
Toep$AcfInput(acf1 + acf2)

for(ii in 1:nSim){
  sim <- Toep$Solve(dX - mean - e1[, ii] - e2[, ii])
  Toep$AcfInput(acf1)
  path.sim[, ii] <- mean + e1[, ii] + Toep$Mult(sim)
}

path.sim.cum <- apply(path.sim, 2, cumsum)

band.sim <- conf.band(path.sim.cum)
rng <- range(c(band.sim, drift))
plot(drift, col = "red", type = "l", ylim = rng, xlab = "", ylab = "")
par(new = TRUE)
plot(band.sim[,1], col = "blue", type = "l", ylim = rng, xlab = "", ylab = "")
par(new = TRUE)
plot(band.sim[,2], col = "blue", type = "l", ylim = rng, xlab = "", ylab = "")


# Drift simulation technique with estimated parameter
dT <- 1/60
N <- 2000
param <- theta.nlm$estimate
mu.para <- param[1]
H.para <- param[2]
lambda.para <- param[3]
Toep <- new(Toeplitz, N)
mean <- mean.fun(mu.para, dT, N)
acf1 <- fbm.acf(H.para, dT, N)
acf2 <- exp2.acf(lambda.para, dT, N)

# simu
nSim <- 100
path.sim <- matrix(NA, N, nSim)
e1 <- rSnorm(nSim, acf = acf1)
e2 <- rSnorm(nSim, acf = acf2)
Toep$AcfInput(acf1 + acf2)

for(ii in 1:nSim){
  sim <- Toep$Solve(dX - mean - e1[, ii] - e2[, ii])
  Toep$AcfInput(acf1)
  path.sim[, ii] <- mean + e1[, ii] + Toep$Mult(sim)
}

path.sim.cum <- apply(path.sim, 2, cumsum)

band.sim <- conf.band(path.sim.cum)
rng <- range(c(band.sim, drift))
plot(drift, col = "red", type = "l", ylim = rng, xlab = "", ylab = "")
par(new = TRUE)
plot(band.sim[,1], col = "blue", type = "l", ylim = rng, xlab = "", ylab = "")
par(new = TRUE)
plot(band.sim[,2], col = "blue", type = "l", ylim = rng, xlab = "", ylab = "")
