# Vignettes Code

# Instruction On FFTW -----------------------------------------------------



# Example of Simulation and Estimation ------------------------------------

setwd("D:/GitHub/SuperGauss/vignettes")
require(numDeriv)
require(SuperGauss)
require(LMN)
require(multiplot)

# Data Simulation ---------------------------------------------------------

# model: 1-d fBM time series
N <- 2000
dT <- 1/60
alpha <- 0.8
acf <- fbm.acf(alpha = alpha, dT = dT, N = N)
dX <- rSnorm(n = 1, acf = acf)
Xt <- cumsum(dX)
tSeq <- 1:N * dT
plot(x = tSeq, y = Xt, col = "blue", type = "l", xlab = "", ylab = "", 
     main = "simulation of 1-d fBM time series")


# Inference, one parameter case, ------------------------------------------
# In order to keep things simple, mu = 0 and Sigma = Identity matrix
# fBM process with parameter alpha: Xt ~ N(0, V(alpha))

# functions for computing the gradiant and hessian of fbm.acf
fbm.acf.grad <- function(alpha, dT, N) {
  if(N == 1) {
    dacf <- dT^alpha * log(alpha)
  } else {
    dacf <- c(0, (dT*(1:N))^alpha * log(dT*(1:N)))
    dacf <- .5 * (dacf[1:N+1] + c(dacf[2], dacf[1:(N-1)]) - 2*dacf[1:N])
  }
  dacf
}

fbm.acf.hess <- function(alpha, dT, N) {
  if(N == 1) {
    d2acf <- dT^alpha * log(alpha)^2
  } else {
    d2acf <- c(0, (dT*(1:N))^alpha * log(dT*(1:N))^2)
    d2acf <- .5 * (d2acf[1:N+1] + c(d2acf[2], d2acf[1:(N-1)]) - 2*d2acf[1:N])
  }
  d2acf
}

# Newton-Raphson Using Function 'nlm'

fbm.loglik <- function(alpha, Xt, dT, acf){
  if(alpha <= 0 || alpha >= 2){
    density <- Inf
  }else{
    dX <- c(Xt[1], diff(Xt))
    N <- length(dX)
    
    acf1 <- fbm.acf(alpha = alpha, dT = dT, N = N)
    acf$setAcf(acf1)
    dacf <- fbm.acf.grad(alpha = alpha, dT = dT, N = N)
    d2acf <- fbm.acf.hess(alpha = alpha, dT = dT, N = N)

    density <- -1 * dSnorm(X = dX, mean = 0, acf = acf, log = TRUE)
    attr(density, "gradient") <- -1 * Snorm.grad(X = dX, mean = 0, acf = acf, dmean = 0, dacf = dacf)
    acf$setAcf(acf1)
    attr(density, "hessian") <- -1 * Snorm.Hess(X = dX, mean = 0, acf = acf, dmean = 0, dacf = dacf, 
                                                d2mean = 0, d2acf = d2acf)
  }
  density
}

acf <- Toeplitz(N)
theta.nlm <- nlm(f = fbm.loglik, p = 1, X = Xt, dT = dT, acf = acf)
signif(c(alpha, theta.nlm$estimate), digits = 3)

# Inference, two parameter case, ------------------------------------------
# In order to keep things simple, mu = 0 and Sigma = Identity matrix
# static fbm error process with parameter alpha and sigma: 
# Xt = fBM(alpha) + sigma * eps, where eps is standard white noise 

# simulating the matern process
alpha <- 0.8
sigma <- 0.2
dT <- 1/30
N <- 2000
acf <- fbm.acf(alpha = alpha, dT = dT, N = N)
dX <- rSnorm(n = 1, acf = acf)
Xt <- cumsum(dX)
Xt <- Xt + sigma * rnorm(N)
tSeq <- 1:N * dT
plot(x = tSeq, y = Xt, col = "blue", type = "l", xlab = "", ylab = "", 
     main = "simulation of 1-d matern time series")

# functions for computing the gradiant and hessian of fbm.acf
static.acf <- function(theta, dT, N){
  alpha <- theta[1]
  sigma <- theta[2]
  acf <- fbm.acf(alpha = alpha, dT = dT, N = N) + sigma^2 * c(2, -1, rep(0, N - 2))
  acf
}

static.acf.grad <- function(theta, dT, N) {
  alpha <- theta[1]
  sigma <- theta[2]
  dacf <- matrix(NA, N, 2)
  if(N == 1) {
    dacf[, 1] <- dT^alpha * log(alpha)
    dacf[, 2] <- 4 * sigma
  } else {
    dacf1 <- c(0, (dT*(1:N))^alpha * log(dT*(1:N)))
    dacf[, 1] <- .5 * (dacf1[1:N+1] + c(dacf1[2], dacf1[1:(N-1)]) - 2*dacf1[1:N])
    dacf[, 2] <- 2 * sigma * c(2, -1, rep(0, N-2))
  }
  dacf
}

static.acf.hess <- function(theta, dT, N) {
  alpha <- theta[1]
  sigma <- theta[2]
  d2acf <- array(0, c(N, 2, 2))
  if(N == 1) {
    d2acf[, 1, 1] <- dT^alpha * log(alpha)^2
    d2acf[, 2, 1] <- 0
    d2acf[, 1, 2] <- 0
    d2acf[, 2, 2] <- 4
  } else {
    d2acf1 <- c(0, (dT*(1:N))^alpha * log(dT*(1:N))^2)
    d2acf[, 1, 1] <- .5 * (d2acf1[1:N+1] + c(d2acf1[2], d2acf1[1:(N-1)]) - 2*d2acf1[1:N])
    d2acf[, 1, 2] <- rep(0, N)
    d2acf[, 2, 1] <- rep(0, N)
    d2acf[, 2, 2] <- 2 * c(2, -1, rep(0, N-2))
  }
  d2acf
}

# Newton-Raphson Using Function 'nlm'

static.loglik <- function(theta, Xt, dT, acf){
  if(theta[1] <= 0 || theta[1] >= 2){
    density <- Inf
  }else{
    dX <- c(Xt[1], diff(Xt))
    N <- length(dX)
    
    acf1 <- static.acf(theta = theta, dT = dT, N = N)
    acf$setAcf(acf1)
    dacf <- static.acf.grad(theta = theta, dT = dT, N = N)
    d2acf <- static.acf.hess(theta = theta, dT = dT, N = N)
    density <- -1 * dSnorm(X = dX, mean = 0, acf = acf, log = TRUE)
    attr(density, "gradient") <- -1 * Snorm.grad(X = dX, mean = 0, acf = acf, dmean = 0, dacf = dacf)
    acf$setAcf(acf1)
    attr(density, "hessian") <- -1 * Snorm.Hess(X = dX, mean = 0, acf = acf, dmean = 0, dacf = dacf, 
                                                d2mean = 0, d2acf = d2acf)
  }
  density
}

static.loglik.alt <- function(theta, Xt, dT, acf){
  if(theta[1] <= 0 || theta[1] >= 2 || theta[2] < 0){
    density <- Inf
  }else{
    dX <- c(Xt[1], diff(Xt))
    N <- length(dX)
    
    acf1 <- static.acf(theta = theta, dT = dT, N = N)
    acf$setAcf(acf1)
    density <- -1 * dSnorm(X = dX, mean = 0, acf = acf, log = TRUE)
  }
  density
}

acf <- Toeplitz(N)
system.time({
  theta.nlm <- nlm(f = static.loglik, p = c(1, 0.1), X = Xt, dT = dT, acf = acf)
})
# time comparison
system.time({
  theta.nlm.alt <- nlm(f = static.loglik.alt, p = c(1, 0.1), X = Xt, dT = dT, acf = acf)
})
theta.nlm$estimate
theta.nlm.alt$estimate
signif(c(alpha, sigma, theta.nlm$estimate), digits = 3)
