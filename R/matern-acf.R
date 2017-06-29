#' @title Autocorrelation(ACF) of Matern processes
#'
#' @description Return the ACF of Matern covariance.
#' @param tseq  length \code{N} vector of timepoints
#' @param lambda timescale
#' @param nu smoothness parameter
#' @return Length \code{N} vector of matern ACF.
#' @details
#' the Matern ACF is given by
#' \deqn{
#' \frac{2^{1-\nu}}{\Gamma(\nu)} \left(\sqrt{2\nu}\frac{t}{\lambda}\right)^\nu K_\nu\left(\sqrt{2\nu} \frac{t}{\lambda}\right)
#' }{
#' 2^(1-v) / \Gamma(v) * (\sqrt{2v} t / \lambda)^v * K_v(\sqrt{2v} t / \lambda)
#' }
#' where K is the Bessel function of second kind.
#' @return An ACF vector of length \code{N}.
#' @examples
#' matern.acf(tseq = 1:10, lambda = 1, nu = 3/2)
#' @export
matern.acf <- function(tseq, lambda, nu) {
  tt <- sqrt(2*nu) * tseq /lambda
  gam <- nu * log(.5 * tt) - lgamma(nu)
  gam <- 2 * exp(gam) * besselK(tt, nu)
  gam[tt == 0] <- 1
  gam
}