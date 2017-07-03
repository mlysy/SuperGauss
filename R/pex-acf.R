#' @title Power-Exponential autocorrelation function
#'
#' @param tseq Length-\code{N} vector of timepoints.
#' @param lambda Timescale parameter.
#' @param rho Power parameter.
#' @return Length-\code{N} autocorrelation vector.
#' @details
#' The power-exponential autocorrelation function is given by:
#' \deqn{
#' \textrm{acf}(t) = \exp \left\{-(\frac{t}{\lambda})^\rho\right\}
#' }{
#' \textrm{acf}(t) = exp (-(t / \lambda)^\rho)
#' }
#' @examples
#' pex.acf(tseq = 1:10, lambda = 1, rho = 2)
#' @export
pex.acf <- function(tseq, lambda, rho) {
  exp(-(tseq/lambda)^rho)
}
