#' @title Convert mean square displacement to autocorrelations.
#'
#' @description Converts the mean squared displacement (MSD) of positions to the autocorrelation (ACF) of the corresponding increments.
#' @param msd Length-\code{N} vector of MSDs at regular timepoints \code{ dt, 2*dt, ..., N*dt}.
#' @return Length \code{N} vector of ACFs.
#' @details
#' For a stationary increments process \eqn{X_t}, converts a sequence \eqn{\eta_1, \ldots, \eta_N} of regularly spaced MSDs,
#' \deqn{
#' \eta_i = E[(X_{i\Delta t} - X_0)^2],
#' }{
#' \eta_i = E[(X_(i*\Delta t) - X_0)^2],
#' }
#' into \eqn{\gamma_1, \ldots, \gamma_N}, a sequence of regularly spaced ACFs,
#' \deqn{
#' \gamma_i = \mathrm{cov}\{X_{(i+1)\Delta t} - X_{i \Delta_i}, X_{\Delta t} - X_{0}\}.
#' }{
#' \gamma_i = cov{X_((i+1)*\Delta t) - X_(i * \Delta t), X_(\Delta t) - X_0}.
#' }
#' This only produces correct results when \code{msd} corresponds to equally-spaced observations.
#' @examples
#' # autocorrelation of fBM increments
#' msd2acf(msd = fbm.msd(tseq = 0:10, H = .3))
#' @export
msd2acf <- function(msd) {
  n <- length(msd)
  gam <- rep(NA, n)
  Gam <- diff(c(0, msd))
  gam[-1] <- .5 * diff(Gam)
  gam[1] <- msd[1]
  gam
}
