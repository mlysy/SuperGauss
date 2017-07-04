#' Matern autocorrelation function.
#'
#' @param tseq Vector of time points at which the autocorrelation is to be calculated.
#' @param lambda Timescale parameter.
#' @param nu Smoothness parameter.
#' @details
#' The Matern autocorrelation is given by
#' \deqn{
#' \mathrm{\scriptsize ACF}(t) = \frac{2^{1-\nu}}{\Gamma(\nu)} \left(\sqrt{2\nu}\frac{t}{\lambda}\right)^\nu K_\nu\left(\sqrt{2\nu} \frac{t}{\lambda}\right),
#' }{
#' acf(t) = 2^(1-\nu)/\Gamma(\nu) * (\sqrt{2\nu} * t/\lambda)^\nu * K_\nu(\sqrt{2\nu} * t/\lambda),
#' }
#' where \eqn{K_\nu(x)} is the modified Bessel function of second kind.
#' @return An autocorrelation vector of length \code{N}.
#' @examples
#' matern.acf(tseq = 1:10, lambda = 1, nu = 3/2)
#' @export
matern.acf <- function(tseq, lambda, nu) {
  tt <- sqrt(2*nu) * abs(tseq)/lambda
  gam <- nu * log(.5 * tt) - lgamma(nu)
  gam <- 2 * exp(gam) * besselK(tt, nu)
  gam[tt == 0] <- 1
  gam
}
