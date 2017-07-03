#--- acf and msd functions -------------------------------------------------

#' Autocorrelation of fBM Increments
#'
#' @param alpha Subdiffusion exponent.
#' @param dT Interobservation time.
#' @param N Number of increment observations.
#' @return An autocorrelation vector of length \eqn{N}.
#' @details
#' The fBM increment autocorrelation is given by:
#' \deqn{\frac{1}{2} \Delta t^\alpha \left[|n-1|^\alpha + |n+1|^\alpha - 2n^\alpha \right]}
#' @examples
#' fbm.acf(alpha = 0.5, dT = 1/60, N = 10)
#' @export
fbm.acf <- function(alpha, dT, N) {
  if(N == 1) {
    acf <- dT^alpha
  } else {
    acf <- (dT*(0:N))^alpha
    acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
  }
  acf
}

#' Autocorrelation of Squared-Exponential
#'
#' @param lambda timescale.
#' @param dT interobservation time.
#' @param N Number of increment observations.
#' @param incr logical; whether or not to return increments.
#' @return An autocorrelation vector of length \eqn{N}.
#' @details
#' The Squared-Exponential autocorrelation is given by:
#' \deqn{\exp \{-(\frac{n\Delta t}{\lambda})^2\}}
#' @examples
#' exp2.acf(lambda = 1, dT = 1/60, N = 200, incr = FALSE)
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
#'
#' @param lambda timescale.
#' @param dT interobservation time.
#' @param N Number of increment observations.
#' @param incr logical; whether or not to return increments.
#' @return An autocorrelation vector of length \eqn{N}.
#' @details
#' The Exponential Autocorrelation is given by:
#' \deqn{\exp \{-\frac{n\Delta t}{\lambda}\}}
#' @examples
#' exp1.acf(lambda = 1, dT = 1/60, N = 200, incr = FALSE)
#' @export
exp1.acf <- function(lambda, dT, N, incr = TRUE) {
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

#' Convert the Position ACF to Increment ACF
#'
#' Converts the autocorrelation of a stationary sequence \eqn{\{X_0, X_1, ..., X_N\}} to those of
#' its increments \eqn{\{X_1-X_0, X_2-X_1, ..., X_N - X_{N-1}\} }.
#' @param gam An autocorrelation sequence of length \eqn{N}.
#' @return An increment autocorrelation sequence of length \eqn{N-1}.
#' @examples
#' acf1 <- runif(10)
#' acf2incr(acf1)
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
#' Converts the MSD of a regularly sampled stationary-increments process \eqn{ \{X_0, X_1, ..., X_N\} }
#' into the ACF of its increments, \eqn{ \{X_1-X_0, X_2-X_1, ..., X_N - X_{N-1}\} }.
#' @param eta the MSD at \eqn{N} regular time points, \eqn{ \{\Delta t, 2\Delta t, ..., N\Delta t\} }.
#' @return the ACF at lags \eqn{ \{0, 1, ..., N-1\} }.
#' @examples
#' msd1 <- runif(10)
#' msd2acf(msd1)
#' @export
msd2acf <- function(eta) {
  N <- length(eta)
  gam <- rep(NA, N)
  Gam <- diff(c(0, eta))
  gam[-1] <- .5 * diff(Gam)
  gam[1] <- eta[1]
  gam
}

#' Convert the Increment ACF to Position MSD
#'
#' Converts the ACF of the increment of a regularly sampled stationary-increments process
#' \eqn{ \{X_1-X_0, X_2-X_1, ..., X_N - X_{N-1} \} } into the MSD of original process \eqn{ \{X_0, X_1, ..., X_N\} }.
#' @param gam the ACF at \eqn{N} lags, \eqn{ \{0, 1, ..., N-1\} }.
#' @return the MSD at regular time points \eqn{ \{\Delta t, 2\Delta t, ..., N\Delta t\} }.
#' @examples
#' acf1 <- runif(10)
#' acf2msd(acf1)
#' @export
acf2msd <- function(gam){
  N <- length(gam)
  eta <- rep(NA, N)
  eta[1] <- gam[1]
  eta[2] <- 2 * (gam[2] + eta[1])
  for(ii in 3:N){
    eta[ii] <- 2 * (gam[ii] + eta[ii - 1]) - eta[ii - 2]
  }
  eta
}

#' MSD of fBM + Dynamic Error
#'
#' @param alpha Subdiffusion exponent
#' @param tau Width of averaging time-window.
#' @param t Vector of time points \eqn{ \{\Delta t, 2\Delta t, ..., N\Delta t\} }
#' @details
#' this function returns the MSD of \eqn{Y_t}, the integral of fBM process \eqn{X_t} with subdiffusion
#' exponent \eqn{\alpha} \deqn{Y_t = \int_{0}^{\tau} X(t-s)ds}. The expression of the MSD is
#' \deqn{\frac{(t+\tau)^\alpha + (t-\tau)^\alpha - 2t^\alpha - 2\tau^\alpha}{(\alpha+1)(\alpha+2)}}
#' @examples
#' fdyn.msd(alpha = 0.8, tau = 1/600, t = (1:200) * 1/60)
#' @export
fdyn.msd <- function(alpha, tau, t){
  tau <- t/tau
  alpha2 <- alpha+2
  eta <- ((tau+1)^alpha2 + (tau-1)^alpha2 - 2*tau^alpha2 - 2)/alpha2
  eta * tau^alpha/(alpha+1)
}

#' ACF of fBM + Dynamic Error Increments
#'
#' @param alpha Subdiffusion exponent
#' @param tau Width of averaging time-window.
#' @param dT interobservation time.
#' @param N Number of increment observations.
#' @details this function returns the autocorrelation of the increment of \eqn{Y_t}, the integral of
#' fBM process \eqn{X_t} with subdiffusion exponent \eqn{\alpha} \deqn{Y_t = \int_{0}^{\tau} X(t-s)ds}
#' @examples
#' fdyn.acf(alpha = 0.8, tau = 1/600, dT = 1/60, N = 200)
#' @export
fdyn.acf <- function(alpha, tau, dT, N) {
  eta <- fdyn.msd(alpha, tau, dT*1:N)
  msd2acf(eta)
}
