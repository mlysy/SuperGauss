#' @title Position mean square displacement(MSD) to increment autocorrelation(ACF)
#'
#' @description Converts the MSD of a regularly sampled stationary-increments process \code{{X0, X1, ..., XN}}
#' into the ACF of its increments, \code{{dX1, dX2, ..., dXN}}
#' @param eta length \code{N} vector of MSD at regular timepoints \code{{dt, 2*dt, ..., N*dt}}
#' @return Length \code{N} vector of ACF.
#' @details 
#' This function only works for evenly-spaced MSD
#' @examples
#' msd2acf(rnorm(10))
#' @export
msd2acf <- function(eta) {
  N <- length(eta)
  gam <- rep(NA, N)
  Gam <- diff(c(0, eta))
  gam[-1] <- .5 * diff(Gam)
  gam[1] <- eta[1]
  gam
}