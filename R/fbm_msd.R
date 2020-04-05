#' Mean square displacement of fractional Brownian motion.
#'
#' @param tseq Length-`N` vector of timepoints.
#' @param H Hurst parameter (between 0 and 1).
#' @return Length-`N` vector of mean square displacements.
#' @details The mean squared displacement (MSD) of a stochastic process `X_t` is defined as
#' ```
#' MSD(t) = E[(X_t - X_0)^2].
#' ```
#' Fractional Brownian motion (fBM) is a continuous Gaussian process with stationary increments, such that its covariance function is entirely defined the MSD, which in this case is `MSD(t) = |t|^(2H)`.
#' @example examples/fbm_msd.R
#' @export
fbm_msd <- function(tseq, H) {
  abs(tseq)^(2*H)
}
