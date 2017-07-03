#' @title Increment autocorrelation(ACF) to Position mean square displacement(MSD).
#'
#' @description Converts the ACF of the increment of a regularly sampled stationary-increments process
#' \code{{dX1, dX2, ..., dXN}} into the MSD of original process \code{{X1, ..., XN}}.
#' @param gam length \code{N} vector of ACF
#' @return Length \code{N} vector of MSD.
#' @examples
#' acf2msd(rnorm(10))
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