# Vignettes Code

# Instruction On FFTW -----------------------------------------------------



# Example of Simulation and Estimation ------------------------------------

setwd("D:/GitHub/SuperGauss/vignettes")
require(numDeriv)
require(SuperGauss)

# Data Simulation ---------------------------------------------------------

# model: 1-d fBM time series
N <- 2000 # length of time series
tseq <- 1:N/60
alpha <- 0.8 # fBM parameter
fbm.acf <- function(tseq, alpha){
  # this function generates the ACF of fBM increments
  msd2acf(fbm.msd(tseq, alpha))
}
acf <- fbm.acf(tseq = tseq, alpha = alpha) # obtain the autocorrelation function of the fBM increments
dX <- rSnorm(n = 2, acf = acf) # simulating the fBM increments
Xt <- apply(dX, 2, cumsum) # recovering the fBM process
plot(x = tseq, y = Xt[, 1], ylim = range(Xt), col = "blue", type = "l", xlab = "", ylab = "", 
     main = "simulation of 2-d fBM time series")
par(new = TRUE)
plot(x = tseq, y = Xt[, 2], ylim = range(Xt), col = "red", type = "l", xlab = "", ylab = "")

# Inference, one parameter case, ------------------------------------------
# In order to keep things simple, mu = 0 and Sigma = Identity matrix
# fBM process with parameter alpha: Xt ~ N(0, V(alpha))

# functions for computing the derivative of the acf of fBM increment
fbm.acf.grad <- function(tseq, alpha) {
  dacf <- c(0, tseq^alpha * log(tseq))
  dacf <- .5 * (dacf[1:N+1] + c(dacf[2], dacf[1:(N-1)]) - 2*dacf[1:N])
  dacf
}

# functions for computing the second derivative of the acf of fBM increment
fbm.acf.hess <- function(tseq, alpha) {
  d2acf <- c(0, tseq^alpha * log(tseq)^2)
  d2acf <- .5 * (d2acf[1:N+1] + c(d2acf[2], d2acf[1:(N-1)]) - 2*d2acf[1:N])
  d2acf
}

# Newton-Raphson Using Function 'nlm'
fbm.loglik <- function(alpha, tseq, Xt, acf){
  if(alpha <= 0 || alpha >= 2){
    density <- Inf
  }else{
    dX <- rbind(Xt[1, ], apply(Xt, 2, diff))
    N <- nrow(dX)
    
    acf1 <- fbm.acf(tseq = tseq, alpha = alpha)
    acf$setAcf(acf1)
    dacf <- fbm.acf.grad(tseq = tseq, alpha = alpha)
    d2acf <- fbm.acf.hess(tseq = tseq, alpha = alpha)
    
    density <- -1 * dSnorm(X = dX, mu = 0, acf = acf, log = TRUE)
    acf$setAcf(acf1)
    attr(density, "gradient") <- -1 * Snorm.grad(X = dX, acf = acf, dacf = dacf)
    acf$setAcf(acf1)
    attr(density, "hessian") <- -1 * Snorm.Hess(X = dX, acf = acf, dacf = dacf, d2acf = d2acf)
  }
  density
}

acf <- Toeplitz(N)
theta.nlm <- nlm(f = fbm.loglik, p = 1, tseq = tseq, X = Xt, acf = acf)$estimate
signif(c(alpha, theta.nlm), digits = 3)

# confidence interval for the MLE
acf1 <- fbm.acf(tseq = tseq, alpha = theta.nlm)
acf$setAcf(acf1)
dacf <- fbm.acf.grad(tseq = tseq, alpha = theta.nlm)
d2acf <- fbm.acf.hess(tseq = tseq, alpha = theta.nlm)
ci.var <- -1 * solve(Snorm.Hess(X = dX, acf = acf, dacf = dacf, d2acf = d2acf))
signif(c(theta.nlm - 1.96 * sqrt(ci.var), theta.nlm + 1.96 * sqrt(ci.var)))

# Inference, two parameter case, ------------------------------------------
# In order to keep things simple, mu = 0 and Sigma = Identity matrix
# static fbm error process with parameter alpha and sigma: 
# Xt = fBM(alpha) + sigma * eps, where eps is standard white noise 

# simulating the matern process
alpha <- 0.8
rho <- 0.2
N <- 2000
tseq <- 1:N/60
acf <- fbm.acf(tseq = tseq, alpha = alpha)
dX <- rSnorm(n = 2, acf = acf)
Xt <- apply(dX, 2, cumsum)
Xt <- Xt + rho * matrix(rnorm(2*N), N, 2)
plot(x = tseq, y = Xt[, 1], ylim = range(Xt), col = "blue", type = "l", xlab = "", ylab = "", 
     main = "simulation of 2-d static error time series")
par(new = TRUE)
plot(x = tseq, y = Xt[, 2], ylim = range(Xt), col = "red", type = "l", xlab = "", ylab = "")


# functions for computing the acf of static error time series increment
static.acf <- function(tseq, theta){
  alpha <- theta[1]
  rho <- theta[2]
  N <- length(tseq)
  acf <- fbm.acf(tseq = tseq, alpha = alpha) + rho^2 * c(2, -1, rep(0, N - 2))
  acf
}

# functions for computing the derivative of the acf of static error time series increment
static.acf.grad <- function(tseq, theta) {
  alpha <- theta[1]
  rho <- theta[2]
  N <- length(tseq)
  dacf <- matrix(NA, N, 2)
  dacf1 <- c(0, tseq^alpha * log(tseq))
  dacf[, 1] <- .5 * (dacf1[1:N+1] + c(dacf1[2], dacf1[1:(N-1)]) - 2*dacf1[1:N])
  dacf[, 2] <- 2 * rho * c(2, -1, rep(0, N-2))
  dacf
}

# functions for computing the secode derivative of the acf of static error time series increment
static.acf.hess <- function(tseq, theta) {
  alpha <- theta[1]
  rho <- theta[2]
  N <- length(tseq)
  d2acf <- array(0, c(N, 2, 2))
  d2acf1 <- c(0, tseq^alpha * log(tseq)^2)
  d2acf[, 1, 1] <- .5 * (d2acf1[1:N+1] + c(d2acf1[2], d2acf1[1:(N-1)]) - 2*d2acf1[1:N])
  d2acf[, 2, 2] <- 2 * c(2, -1, rep(0, N-2))
  d2acf
}

# Newton-Raphson Using Function 'nlm'

static.loglik <- function(theta, tseq, Xt, acf){
  if(theta[1] <= 0 || theta[1] >= 2){
    density <- Inf
  }else{
    dX <- rbind(Xt[1, ], apply(Xt, 2, diff))
    N <- nrow(dX)
    
    acf1 <- static.acf(tseq = tseq, theta = theta)
    dacf <- static.acf.grad(tseq = tseq, theta = theta)
    d2acf <- static.acf.hess(tseq = tseq, theta = theta)
    acf$setAcf(acf1)
    density <- -1 * dSnorm(X = dX, acf = acf, log = TRUE)
    attr(density, "gradient") <- -1 * Snorm.grad(X = dX, acf = acf, dacf = dacf)
    acf$setAcf(acf1)
    attr(density, "hessian") <- -1 * Snorm.Hess(X = dX, acf = acf, dacf = dacf, d2acf = d2acf)
  }
  density
}

acf <- Toeplitz(N)
theta.nlm <- nlm(f = static.loglik, p = c(1, 0.1), tseq = tseq, X = Xt, acf = acf)$estimate
signif(c(alpha, rho, theta.nlm), digits = 3)

# confidence interval for the MLE
acf1 <- static.acf(tseq = tseq, theta = theta.nlm)
dacf <- static.acf.grad(tseq = tseq, theta = theta.nlm)
d2acf <- static.acf.hess(tseq = tseq, theta = theta.nlm)
acf$setAcf(acf1)
ci.var <- -1 * solve(Snorm.Hess(X = dX, acf = acf, dacf = dacf, d2acf = d2acf))
signif(c(theta.nlm[1] - 1.96 * sqrt(ci.var[1, 1]), theta.nlm[1] + 1.96 * sqrt(ci.var[1, 1])))
signif(c(theta.nlm[2] - 1.96 * sqrt(ci.var[2, 2]), theta.nlm[2] + 1.96 * sqrt(ci.var[2, 2])))
