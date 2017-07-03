#' @title Autocorrelation(ACF) of power exponential processes
#'
#' @description Return the ACF of power exponential process at different lags
#' @param tseq length \code{N} vector of timepoints
#' @param lambda Timescale parameter
#' @param rho Power parameter
#' @return Length \code{N} vector of power-exponential ACF.
#' @details
#' The power-exponential ACF is given by:
#' \deqn{
#' \textrm{acf}(t) = \exp \{-(\frac{t}{\lambda})^\rho\}
#' }{
#' \textrm{acf}(t) = exp (-(t / \lambda)^\rho)
#' }
#' @examples
#' pex.acf(tseq = 1:10, lambda = 1, rho = 2)
#' @export
pex.acf <- function(tseq, lambda, rho) {
  exp(-(tseq/lambda)^rho)
}