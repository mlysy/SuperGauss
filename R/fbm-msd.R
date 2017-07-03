#' Mean square displacement of fractional Brownian motion.
#'
#' @param tseq Length-\code{N} vector of timepoints.
#' @param H Hurst parameter (between 0 and 1).
#' @return Length \code{N} vector of mean square displacements.
#' @details The mean squared displacement (MSD) of a stochastic process $X_t$ is defined as
#' \deqn{
#' \mathrm{\scriptsize MSD}_X(t) = E[(X_t - X(0)_0)^2].
#' }{
#' MSD_X(t) = E[(X_t - X_0)^2].
#' }
#' Fractional Brownian motion (fBM) is a continuous Gaussian process with stationary increments, such that its covariance function is entirely defined the MSD, which in this case is \eqn{\textrm{\small MSD}_X(t) = |t|^{2H}}{MSD_X(t) = |t|^(2H)}.
#' @examples
#' fbm.msd(tseq = 1:10, H = 0.4)
#' @export
fbm.msd <- function(tseq, H) {
  abs(tseq)^(2*H)
}
