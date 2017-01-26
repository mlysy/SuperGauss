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
d <- 1
dT <- 1/60
H <- 0.4
lambda <- 200 * dT
mu <- 1
theta <- c(mu, H, lambda)

# three parts of X
drift <- mean.fun(mu, dT, N)
fbm <- rSnorm(d, acf = fbm.acf(H, dT, N))
exp2 <- rSnorm(d, acf = exp2.acf(lambda, dT, N))

# X and plot
dX <- matrix(drift + fbm + exp2, N)
X <- cumsum(dX)
plot(X, col = "blue")

# Estimating Parameters
# using Newton-Raphson

# the acf function and its derivatives
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

# dSnorm for checking
dSnorm.check <- function(theta, X, dT, Toep){
  n <- nrow(X)
  mu <- theta[1]
  H <- theta[2]
  lambda <- theta[3]
  mean <- mean.fun(mu, dT, n)
  acf <- acf.fun(H, lambda, dT, n)
  X <- X - mean
  Toep$AcfInput(acf)
  density <- crossprod(X, Toep$Solve(X))
  density <- density + Toep$Det()
  density / -2
}


# parameter theta = {H, lambda}
Newton.Raphson <- function(theta, X, dT, Toep, check = FALSE){
  n <- nrow(X)
  d <- ncol(X)
  p <- 3
  mu <- theta[1]
  H <- theta[2]
  lambda <- theta[3]
  mean <- mean.fun(mu, dT, n)
  acf <- acf.fun(H, lambda, dT, N)
  
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

  grad <- Snorm.grad(X, mean, acf, dmean, dacf, Toep)
  hess <- Snorm.Hess(X, mean, acf, dmean, dacf, d2mean, d2acf, Toep)
  
  # checking part
  if(check){
    grad.check <- grad(func = dSnorm.check, x = theta, X = X, dT = dT, Toep = Toep)
    hess.check <- hessian(func = dSnorm.check, x = theta, X = X, dT = dT, Toep = Toep)
    
    if(max(range(c(grad - grad.check, hess - hess.check))) > 1e-4){
      print(hess)
      print(hess.check) 
    }
  }
  # check end
  
  theta.new <- theta - solve(hess, grad)
  theta.new
}

Newton.Raphson(c(1, 0.5, 3), dX, dT, Toep, TRUE)

iterate <- function(X, dT, start, Toep, check = FALSE, difference = 1e-4){
  theta.old <- start
  theta.new <- Newton.Raphson(theta.old, X, dT, Toep)
  loop <- 1
  while(max(theta.new - theta.old) > difference){
    theta.old <- theta.new
    theta.new <- Newton.Raphson(theta.old, X, dT, Toep, check)
    loop <- loop + 1
  }
  list(est = theta.new, times = loop)
}

Toep <- new(Toeplitz, N)
theta.start <- c(0.7, 0.5, 4)
theta.est <- iterate(dX, dT, theta.start, Toep)
signif(cbind(theta.est$est, theta))
# simulating drift with estimated data

rep <- 100

mu.est <- theta.est$est[1]
H.est <- theta.est$est[2]
lambda.est <- theta.est$est[3]

dX.sim <- matrix(NA, N, rep)
for(ii in 1:rep){
  drift.tmp <- mean.fun(mu.est, dT, N)
  fbm.tmp <- rSnorm(1, acf = fbm.acf(H.est, dT, N))
  exp2.tmp <- rSnorm(1, acf = exp2.acf(lambda.est, dT, N))
  dX.sim[, ii] <- drift.tmp + fbm.tmp + exp2.tmp
}
X.sim <- apply(dX.sim, 2, cumsum)

band <- conf.band(X.sim)

plot(X, col = "red", ylim = c(-2, 46))
par(new = TRUE)
plot(band[,1], col = "green", ylim = c(-2, 46))
par(new = TRUE)
plot(band[,2], col = "green", ylim = c(-2, 46))

## true parameter simulation

dX.true <- matrix(NA, N, rep)
for(ii in 1:rep){
  drift.tmp <- mean.fun(mu, dT, N)
  fbm.tmp <- rSnorm(1, acf = fbm.acf(H, dT, N))
  exp2.tmp <- rSnorm(1, acf = exp2.acf(lambda, dT, N))
  dX.true[, ii] <- drift.tmp + fbm.tmp + exp2.tmp
}
X.true <- apply(dX.true, 2, cumsum)

band.true <- conf.band(X.true)

plot(X, col = "red", ylim = c(-2, 46))
par(new = TRUE)
plot(band[,1], col = "green", ylim = c(-2, 46))
par(new = TRUE)
plot(band[,2], col = "green", ylim = c(-2, 46))