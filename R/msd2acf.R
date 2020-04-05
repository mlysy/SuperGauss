#' Convert mean square displacement of positions to autocorrelation of increments.
#'
#' Converts the mean squared displacement (MSD) of a stationary increments sequence `x = (x_0, x_1, ..., x_N)` positions to the autocorrelation of the corresponding increments `dx = (x_1 - x_0, ..., x_N - x_(N-1))`.
#'
#' @param msd Length-`N` MSD vector, i.e., excluding `x_0` which is assumed to be zero.
#' @return Length-`N` autocorrelation vector.
#' @example examples/msd2acf.R
#' @export
msd2acf <- function(msd) {
  n <- length(msd)
  gam <- rep(NA, n)
  Gam <- diff(c(0, msd))
  gam[-1] <- .5 * diff(Gam)
  gam[1] <- msd[1]
  gam
}
