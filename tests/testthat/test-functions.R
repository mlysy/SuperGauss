test_fbm_acf <- function(alpha, dT, N) {
  if(N == 1) {
    acf <- dT^alpha
  } else {
    acf <- (dT*(0:N))^alpha
    acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
  }
  acf
}

test_exp_acf <- function(lambda, order, N) {
  # process autocorrelation
  gam <- exp(-(0:N/N)^order * lambda)
  ans <- gam[1:N]
  ans
}

test_matern_acf <- function(lambda, nu, N) {
  # process autocorrelation
  tt <- sqrt(2*nu) * (0:N) / lambda
  gam <- nu * log(.5 * tt) - lgamma(nu)
  gam <- 2 * exp(gam) * besselK(tt, nu)
  gam[tt == 0] <- 1
  ans <- gam[1:N]
  ans
}

test_acf_func <- function(N, type, first0 = FALSE) {
  lambda <- 1 / runif(1, 1, 3)
  alpha <- runif(1, .2, .9)
  nu <- 3
  dT <- runif(1, .1, .9)
  type <- match.arg(type,
                    choices = c("exp", "exp2", "fbm", "matern", "zero", "rnd"))
  if(type == "exp2") {
    acf <- test_exp_acf(lambda, 2, N)
  } else if(type == "exp") {
    acf <- test_exp_acf(lambda, 1, N)
  } else if(type == "fbm") {
    acf <- test_fbm_acf(alpha, dT, N)
  } else if(type == "matern") {
    acf <- test_matern_acf(lambda, nu, N)
  } else if(type == "zero") {
    acf <- rep(0, N)
  } else {
    acf <- rnorm(N)
  }
  if(first0) acf[1] <- 0
  acf
}

test_fbm_acf_grad <- function(alpha, dT, N) {
  if(N == 1) {
    dacf <- dT^alpha * log(alpha)
  } else {
    dacf <- c(0, (dT*(1:N))^alpha * log(dT*(1:N)))
    dacf <- .5 * (dacf[1:N+1] + c(dacf[2], dacf[1:(N-1)]) - 2*dacf[1:N])
  }
  dacf
}

test_fbm_acf_hess <- function(alpha, dT, N) {
  if(N == 1) {
    d2acf <- dT^alpha * log(alpha)^2
  } else {
    d2acf <- c(0, (dT*(1:N))^alpha * log(dT*(1:N))^2)
    d2acf <- .5 * (d2acf[1:N+1] + c(d2acf[2], d2acf[1:(N-1)]) - 2*d2acf[1:N])
  }
  d2acf
}


test_drift_func <- function(mu, N){
  rep(mu^2, N)
}

test_drift_grad <- function(mu, N){
  rep(2 * mu, N)
}

test_drift_hess <- function(mu, N){
  rep(2, N)
}
