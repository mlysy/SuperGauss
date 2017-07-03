#' Superfast inference for stationary Gaussian time series.
#'
#' @details While likelihood calculations with stationary Gaussian time series generally scale as \eqn{\mathcal O(N^2)}{O(N^2)} in the number of observations, this package implements an algorithm which scales as \eqn{\mathcal O(N \log^2 N)}{O(N \log^2 N)}.  "Superfast" algorithms for loglikelihood gradients and Hessians are also provided.  The underlying C++ code is distributed through a header-only library found in the installed package's \code{include} directory.
#' @examples
#' # Superfast inference for the timescale parameter of
#' # the exponential autocorrelation function
#' exp.acf <- function(lambda) exp(-(1:N-1)/lambda)
#'
#' # simulate data
#' lambda0 <- 1
#' N <- 1000
#' X <- rSnorm(n = 1, acf = exp.acf(lambda0))
#'
#' # loglikelihood function
#' Toep <- Toeplitz(n = N) # allocate memory for a Toeplitz matrix object
#' loglik <- function(lambda) {
#'   Toep$setAcf(acf = exp.acf(lambda))
#'   dSnorm(X = X, acf = Toep, log = TRUE)
#' }
#'
#' # maximum likelihood estimation
#' optimize(f = loglik, interval = c(.2, 5), maximum = TRUE)
#' @docType package
#' @name SuperGauss
#' @importFrom Rcpp evalCpp
#' @importFrom methods new show
#' @importFrom stats rnorm
#' @useDynLib SuperGauss, .registration = TRUE
NULL
