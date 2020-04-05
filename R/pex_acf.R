#' Power-exponential autocorrelation function.
#'
#' @param tseq Vector of `N` time points at which the autocorrelation is to be calculated.
#' @param lambda Timescale parameter.
#' @param rho Power parameter.
#' @return An autocorrelation vector of length `N`.
#'
#' @details
#' The power-exponential autocorrelation function is given by:
#' \deqn{
#' \mathrm{\scriptsize ACF}(t) = \exp \left\{-(t/\lambda)^\rho\right\}.
#' }{
#' acf(t) = exp (-(t / \lambda)^\rho).
#' }
#' @example examples/pex_acf.R
#' @export
pex_acf <- function(tseq, lambda, rho) {
  exp(-abs(tseq/lambda)^rho)
}
