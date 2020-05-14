library(SuperGauss)

test_fbm_acf <- function(alpha, dT, N) {
  msd <- fbm_msd(dT*(1:N), alpha/2)
  msd2acf(msd)
  ## if(N == 1) {
  ##   acf <- dT^alpha
  ## } else {
  ##   acf <- (dT*(0:N))^alpha
  ##   acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
  ## }
  ## acf
}

## alpha <- runif(1, 0, 2)
## dT <- runif(1, 0, 5)
## N <- 1
## test_fbm_acf(alpha, dT, N) - test_fbm_acf2(alpha, dT, N)

test_exp_acf <- function(lambda, order, N) {
  # process autocorrelation
  ## gam <- exp(-(0:N/N)^order * lambda)
  ## ans <- gam[1:N]
  ## ans
  pex_acf(0:(N-1)/N, rho = order, lambda = 1/lambda^(1/order))
}

## lambda <- 1 / runif(1, 1, 3)
## order <- runif(1, 0, 5)
## N <- 5
## test_exp_acf(lambda, order, N) - test_exp_acf2(lambda, order, N)

test_matern_acf <- function(lambda, nu, N) {
  # process autocorrelation
  ## tt <- sqrt(2*nu) * (0:N) / lambda
  ## gam <- nu * log(.5 * tt) - lgamma(nu)
  ## gam <- 2 * exp(gam) * besselK(tt, nu)
  ## gam[tt == 0] <- 1
  ## ans <- gam[1:N]
  ## ans
  matern_acf(0:(N-1), lambda, nu)
}

## lambda <- 1 / runif(1, 1, 3)
## nu <- runif(1, 0, 5)
## N <- 5
## test_matern_acf(lambda, nu, N) - test_matern_acf2(lambda, nu, N)


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

# log-density of z ~ NormalToeplitz(gamma)
toep_ldens <- function(z, gamma) {
  mvtnorm::dmvnorm(x = z, sigma = toeplitz(gamma), log = TRUE)
}


#--- circulant functions -------------------------------------------------------

# columnwise fft for real inputs
fft <- function(x) {
  eV <- is.vector(x)
  x <- as.matrix(x)
  if(nrow(x) == 1) {
    y <- x + 0i
  }  else {
    y <- apply(as.matrix(x), 2, fftw::FFT)
  }
  if(eV) y <- drop(y)
  y
}

# columnwise normalized inverse fft for real input
ifft <- function(x) {
  eV <- is.vector(x)
  x <- as.matrix(x)
  if(nrow(x) == 1) {
    y <- x + 0i
  }  else {
    y <- apply(as.matrix(x), 2, fftw::FFT, inverse = TRUE)
  }
  y <- Re(y)/nrow(y)
  if(eV) y <- drop(y)
  y
}


# convert uacf to full acf.
unfold_acf <- function(N, uacf) {
  n <- length(uacf)
  if(n != floor(N/2) + 1) stop("uacf has wrong length.")
  acf <- rep(NA, N)
  acf[1:n] <- uacf
  if(N > 1) {
    eN <- (2*n) == (N+2)
    id <- n - eN + (2:n-1)
    acf[n - eN + (2:n-1)] <- uacf[n:2]
  }
  acf
}

#' Circulant matrix with first row `x`.
circulant <- function(x) {
  N <- length(x)
  t(sapply(1:N-1, function(ii) x[(1:N-1 - ii) %% N + 1]))
}

# log-density of z ~ NormalCirculant(nu)
circ_ldens <- function(z, nu) {
  mvtnorm::dmvnorm(z, log = TRUE,
                   sigma = toeplitz(unfold_acf(length(z), nu)))
}
