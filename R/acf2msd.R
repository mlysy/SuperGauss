#' Convert autocorrelation of stationary increments to mean squared displacement of posititions.
#'
#' Converts the autocorrelation (ACF) of stationary increments to the mean squared displacement (MSD) of the corresponding positions.
#' @param acf length-\code{N} ACF vector of a stationary increment sequence.
#' @return Length-\code{N} MSD vector of the corresponding positions.
#' @details If \eqn{X(t)} is a stationary increments process, then \eqn{\Delta X_0, \Delta X_1, \ldots} with
#' \deqn{
#' \Delta X_n = X((n+1)\Delta t) - X(n \Delta t)
#' }
#' is a stationary time series.  This function converts the ACF of this series into the MSD of the corresponding positions, namely returns the sequence \eqn{\eta_1, \ldots, \eta_N}, where \eqn{\eta_i = \mathrm{var}(X(i\Delta t))}{\eta_i = var(X(i\Delta t))}.
#' @examples
#' acf2msd(acf = exp(-(0:10)))
#' @export
acf2msd <- function(acf){
  N <- length(acf)
  msd <- rep(NA, N)
  msd[1] <- acf[1]
  msd[2] <- 2 * (acf[2] + msd[1])
  for(ii in 3:N){
    msd[ii] <- 2 * (acf[ii] + msd[ii - 1]) - msd[ii - 2]
  }
  msd
}
