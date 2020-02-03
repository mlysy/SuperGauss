#' @title Gradient of the loglikelihood of a multivariate normal with Toeplitz variance matrix.
#'
#' @description Superfast evaluation of loglikelihood gradient.
#' @param X A length-\code{N} vector of multivariate normal observations.
#' @param mu A scalar or length-\code{N} vector of means.  If missing defaults to the vector of zeros.
#' @param acf A \code{Toeplitz} object or length-\code{N} vector containing the first column of the Toeplitz variance matrix.
#' @param dmu A length-\code{p} vector or \code{N x p} matrix of partial derivatives of \code{mu} along the columns.  If missing defaults to a matrix of zeros.
#' @param dacf An \code{N x p} matrix with the partial derivatives of \code{acf} along the columns.
#' @return A length-\code{p} vector containing the gradient of the loglikelihood.
#' @examples
#' # two parameter inference
#' acf.fun <- function(theta) theta[2]^2 * exp(-theta[1]*(1:N-1))
#' mu.fun <- function(theta) theta[1]
#'
#' # partial derivatives
#' dacf.fun <- function(theta) {
#'   ea <- exp(-theta[1]*(1:N-1))
#'   cbind(-theta[1]*theta[2]^2 * ea, 2*theta[2] * ea)
#' }
#' dmu.fun <- function(theta) c(1, 0)
#'
#' # generate data
#' N <- 300
#' theta <- rexp(2)
#' X <- rSnorm(n = 1, acf = acf.fun(theta)) + mu.fun(theta)
#'
#' # likelihood gradient
#' Snorm.grad(X = X, mu = mu.fun(theta), dmu = dmu.fun(theta),
#'            acf = acf.fun(theta), dacf = dacf.fun(theta))
#' @export
Snorm.grad <- function(X, mu, acf, dmu, dacf) {
  # format arguments
  if(!is.vector(X)) stop("X must be a vector.")
  N <- length(X)
  acf <- .format.acf(acf, N)
  p <- .get.p(dmu, dacf)
  ## if(!is.matrix(dacf) || nrow(dacf) != N) {
  ##   stop("dacf must be a matrix with nrow(dacf) == length(X).")
  ## }
  ## p <- ncol(dacf)
  Mu <- .format.mu(mu = mu, dmu = dmu, N = N, p = p)
  mu <- Mu$mu
  dmu <- Mu$dmu
  dacf <- .format.dacf(dacf = dacf, N = N, p = p)$dacf
  # gradient calculation
  Z <- c(X - mu)
  SigZ <- c(solve(acf, Z))
  grad <- rep(NA, p)
  for(ii in 1:p) {
    grad[ii] <- crossprod(SigZ, toep.mult(dacf[,ii], SigZ))
    grad[ii] <- grad[ii] - acf$traceT2(dacf[,ii])
  }
  .5 * grad + colSums(dmu * SigZ)
}
