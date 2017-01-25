#--- acf and msd functions -------------------------------------------------

#' Autocorrelation of fBM Increments
#'
#' @param H Hurst parameter.  Note that subdiffusion parameter \code{alpha = 2H}.
#' @param dT Interobservation time.
#' @param N Number of increment observations.
#' @return An autocorrelation vector of length \code{N}.
#' @details The fBM autocorrelation is given by:
#' @export
fbm.acf <- function(H, dT, N, incr = TRUE){
   gam <- (dT*(0:N+1))^(2*H)
   if(incr){
     # increments
     ans <- -1/2 * acf2incr(gam)
   } else{
     # observations
     ans <- gam[1:N]
   }
   ans
}

#' Autocorrelation of Squared-Exponential
#'
#' @param lambda timescale.
#' @param dT interobservation time.
#' @param N Number of increment observations.
#' @param incr logical; whether or not to return increments.
#' @return An autocorrelation vector of length \code{N}.
#' @details The Squared-Exponential autocorrelation is given by
#' @export
exp2.acf <- function(lambda, dT, N, incr = TRUE) {
  # process autocorrelation
  gam <- exp(-(0:N*dT/lambda)^2)
  if(incr) {
    # increments
    ans <- acf2incr(gam)
  } else {
    # observations
    ans <- gam[1:N]
  }
  ans
}

#' Autocorrelation of Exponential
#' @param lambda timescale.
#' @param dT interobservation time.
#' @param N Number of increment observations.
#' @param incr logical; whether or not to return increments.
#' @return An autocorrelation vector of length \code{N}.
#' @details The Exponential Autocorrelation is given by
#' @export
exp.acf <- function(lambda, dT, N, incr = TRUE) {
  # process autocorrelation
  gam <- exp(-(0:N*dT/lambda))
  if(incr) {
    # increments
    ans <- acf2incr(gam)
  } else {
    # observations
    ans <- gam[1:N]
  }
  ans
}

#' Autocorrelation of Matern
#'
#' @param lambda timescale.
#' @param nu smoothness parameter.
#' @param dT interobservation time.
#' @param N Number of increment observations.
#' @details the Matern autocorrelation is
#' @return An autocorrelation vector of length \code{N}.
#' @export
matern.acf <- function(lambda, nu, dT, N, incr = TRUE) {
  # process autocorrelation
  tt <- sqrt(2*nu) * (0:N)*dT/lambda
  gam <- nu * log(.5 * tt) - lgamma(nu)
  gam <- 2 * exp(gam) * besselK(tt, nu)
  gam[tt == 0] <- 1
  if(incr) {
    # increments
    ans <- acf2incr(gam)
  } else {
    # observations
    ans <- gam[1:N]
  }
  ans
}

#' Convert the Position ACF to Increment ACF
#'
#' Converts the autocorrelation of a stationary sequence \code{X0, X1, ..., XN} to those of its increments \code{X1-X0, X2-X1, ..., XN - X(N-1)}.
#' @param gam An autocorrelation sequence of length \code{nObs}.
#' @return An increment autocorrelation sequence of length \code{nObs-1}.
#' @export
acf2incr <- function(gam) {
  N <- length(gam)-1
  if(N == 1) {
    igam <- 2*(gam[1]-gam[2])
  } else {
    igam <- 2*gam[1:N] - gam[1:N+1] - gam[c(2, 1:(N-1))]
  }
  igam
}

#' Convert the Position MSD to Increment ACF
#'
#' Converts the MSD of a regularly sampled stationary-increments process \code{X0, X1, ..., XN} into the ACF of its increments, \code{X1-X0, X2-X1, ..., XN - X(N-1)}.
#' @param eta the MSD at \code{nObs} regular time points, \code{dT, 2*dT, ..., nObs*dT}.
#' @return the ACF at lags \code{0, 1, ..., nObs-1}.
#' @export
msd2acf <- function(eta) {
  N <- length(eta)
  gam <- rep(NA, N)
  Gam <- diff(c(0, eta))
  gam[-1] <- .5 * diff(Gam)
  gam[1] <- eta[1]
  gam
}

#' MSD of fBM + Dynamic Error
#'
#' @param alpha Subdiffusion exponent
#' @param sigma Width of averaging time-window.
#' @param t Vector of time points
#' @export
fdyn.msd <- function(alpha, sigma, t){
  tau <- t/sigma
  alpha2 <- alpha+2
  eta <- ((tau+1)^alpha2 + (tau-1)^alpha2 - 2*tau^alpha2 - 2)/alpha2
  eta * sigma^alpha/(alpha+1)
}

#' ACF of fBM + Dynamic Error Increments
#'
#' @export
fdyn.acf <- function(alpha, sigma, dT, N) {
  eta <- fdyn.msd(alpha, sigma, dT*1:N)
  msd2acf(eta)
}

mean.fun <- function(mu, dT, N){
  mu^2 * dT * matrix(1, N, 1)
}