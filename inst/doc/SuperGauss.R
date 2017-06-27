## ------------------------------------------------------------------------
require(SuperGauss)
N <- 2000 # length of time series
dT <- 1/60 # frequency of observation
alpha <- 0.8 # fBM parameter
acf <- fbm.acf(alpha = alpha, dT = dT, N = N) # obtain the autocorrelation function of the fBM increments
dX <- rSnorm(n = 2, acf = acf) # simulating the fBM increments
Xt <- apply(dX, 2, cumsum) # recovering the fBM process
tSeq <- 1:N * dT
plot(x = tSeq, y = Xt[, 1], ylim = range(Xt), col = "blue", type = "l", xlab = "", ylab = "", 
     main = "simulation of 2-d fBM time series")
par(new = TRUE)
plot(x = tSeq, y = Xt[, 2], ylim = range(Xt), col = "red", type = "l", xlab = "", ylab = "")

## ------------------------------------------------------------------------
# functions for computing the derivative of the acf of fBM increment
fbm.acf.grad <- function(alpha, dT, N) {
  if(N == 1) {
    dacf <- dT^alpha * log(alpha)
  } else {
    dacf <- c(0, (dT*(1:N))^alpha * log(dT*(1:N)))
    dacf <- .5 * (dacf[1:N+1] + c(dacf[2], dacf[1:(N-1)]) - 2*dacf[1:N])
  }
  dacf
}

# functions for computing the second derivative of the acf of fBM increment
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
    dX <- apply(Xt, 2, diff)
    N <- nrow(dX)
    
    acf1 <- fbm.acf(alpha = alpha, dT = dT, N = N)
    acf$setAcf(acf1)
    dacf <- fbm.acf.grad(alpha = alpha, dT = dT, N = N)
    d2acf <- fbm.acf.hess(alpha = alpha, dT = dT, N = N)

    density <- -1 * dSnorm(X = dX, mu = 0, acf = acf, log = TRUE)
    acf$setAcf(acf1)
    attr(density, "gradient") <- -1 * Snorm.grad(X = dX, acf = acf, dacf = dacf)
    acf$setAcf(acf1)
    attr(density, "hessian") <- -1 * Snorm.Hess(X = dX, acf = acf, dacf = dacf, d2acf = d2acf)
  }
  density
}

acf <- Toeplitz(N - 1)
theta.nlm <- nlm(f = fbm.loglik, p = 1, X = Xt, dT = dT, acf = acf)
signif(c(alpha, theta.nlm$estimate), digits = 3)

## ------------------------------------------------------------------------
# simulating the matern process
alpha <- 0.8
rho <- 0.2
dT <- 1/30
N <- 2000
acf <- fbm.acf(alpha = alpha, dT = dT, N = N)
dX <- rSnorm(n = 2, acf = acf)
Xt <- apply(dX, 2, cumsum)
Xt <- Xt + rho * matrix(rnorm(2*N), N, 2)
tSeq <- 1:N * dT
plot(x = tSeq, y = Xt[, 1], ylim = range(Xt), col = "blue", type = "l", xlab = "", ylab = "", 
     main = "simulation of 2-d static error time series")
par(new = TRUE)
plot(x = tSeq, y = Xt[, 2], ylim = range(Xt), col = "red", type = "l", xlab = "", ylab = "")


# functions for computing the acf of static error time series increment
static.acf <- function(theta, dT, N){
  alpha <- theta[1]
  rho <- theta[2]
  acf <- fbm.acf(alpha = alpha, dT = dT, N = N) + rho^2 * c(2, -1, rep(0, N - 2))
  acf
}

# functions for computing the derivative of the acf of static error time series increment
static.acf.grad <- function(theta, dT, N) {
  alpha <- theta[1]
  rho <- theta[2]
  dacf <- matrix(NA, N, 2)
  if(N == 1) {
    dacf[, 1] <- dT^alpha * log(alpha)
    dacf[, 2] <- 4 * rho
  } else {
    dacf1 <- c(0, (dT*(1:N))^alpha * log(dT*(1:N)))
    dacf[, 1] <- .5 * (dacf1[1:N+1] + c(dacf1[2], dacf1[1:(N-1)]) - 2*dacf1[1:N])
    dacf[, 2] <- 2 * rho * c(2, -1, rep(0, N-2))
  }
  dacf
}

# functions for computing the secode derivative of the acf of static error time series increment
static.acf.hess <- function(theta, dT, N) {
  alpha <- theta[1]
  rho <- theta[2]
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
    dX <- apply(Xt, 2, diff)
    N <- nrow(dX)
    
    acf1 <- static.acf(theta = theta, dT = dT, N = N)
    dacf <- static.acf.grad(theta = theta, dT = dT, N = N)
    d2acf <- static.acf.hess(theta = theta, dT = dT, N = N)
    acf$setAcf(acf1)
    density <- -1 * dSnorm(X = dX, acf = acf, log = TRUE)
    attr(density, "gradient") <- -1 * Snorm.grad(X = dX, acf = acf, dacf = dacf)
    acf$setAcf(acf1)
    attr(density, "hessian") <- -1 * Snorm.Hess(X = dX, acf = acf, dacf = dacf, d2acf = d2acf)
  }
  density
}

acf <- Toeplitz(N - 1)
theta.nlm <- nlm(f = static.loglik, p = c(1, 0.1), X = Xt, dT = dT, acf = acf)
signif(c(alpha, rho, theta.nlm$estimate), digits = 3)

