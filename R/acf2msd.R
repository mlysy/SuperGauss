#' Convert autocorrelation of stationary increments to mean squared displacement of posititions.
#'
#' Converts the autocorrelation of a stationary increment sequence `dx = (x_1 - x_0, ..., x_N - x_(N-1))` to the mean squared displacement (MSD) of the corresponding positions, i.e., `MSD_i = E[(x_i - x_0)^2]`.
#'
#' @param acf Length-`N` autocorrelation vector of a stationary increment sequence.
#' @return Length-`N` MSD vector of the corresponding positions.
#' @example examples/acf2msd.R
#'
#' @export
acf2msd <- function(acf) {
  N <- length(acf)
  msd <- rep(NA, N)
  msd[1] <- acf[1]
  msd[2] <- 2 * (acf[2] + msd[1])
  for(ii in 3:N) {
    msd[ii] <- 2 * (acf[ii] + msd[ii - 1]) - msd[ii - 2]
  }
  msd
}
