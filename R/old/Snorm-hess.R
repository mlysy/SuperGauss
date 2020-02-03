#' @title Hessian of the loglikelihood of a multivariate normal with Toeplitz variance matrix.
#' @description Superfast evaluation of loglikelihood Hessian.
#' @inheritParams Snorm.grad
#' @param d2mu A \code{p x p} matrix or \code{N x p x p} array of second partial derivatives of \code{mu}.  If missing defaults to zeros.
#' @param d2acf A \code{N x p x p} array of second partial derivatives of \code{acf}.
#' @return The \code{p x p} Hessian matrix of the loglikelihood.
#' @examples
#' # two parameter inference
#' acf.fun <- function(theta) theta[2]^2 * exp(-(1:N-1))
#' mu.fun <- function(theta) theta[1] * (1:N) + log(theta[2] + 1:N)
#'
#' # partial derivatives
#' dacf.fun <- function(theta) {
#'   cbind(0, 2*theta[2] * exp(-(1:N-1)))
#' }
#' dmu.fun <- function(theta) cbind(1:N, 1/(theta[2] + 1:N))
#'
#' # 2nd order partials
#' d2acf.fun <- function(theta) {
#'   H <- array(0, dim = c(N, 2, 2))
#'   H[,2,2] <- 2*exp(-(1:N-1))
#'   H
#' }
#' d2mu.fun <- function(theta) {
#'   H <- array(0, dim = c(N, 2, 2))
#'   H[,2,2] <- -1/(theta[2] + 1:N)^2
#'   H
#' }
#'
#' # generate data
#' N <- 300
#' theta <- rexp(2)
#' X <- rSnorm(n = 1, acf = acf.fun(theta)) + mu.fun(theta)
#'
#' # likelihood Hessian
#' Snorm.hess(X = X, mu = mu.fun(theta), acf = acf.fun(theta),
#'            dmu = dmu.fun(theta), dacf = dacf.fun(theta),
#'            d2mu = d2mu.fun(theta), d2acf = d2acf.fun(theta))
#' @export
Snorm.hess <- function(X, mu, acf, dmu, dacf, d2mu, d2acf) {
  if(!is.vector(X)) stop("X must be a vector.")
  N <- length(X)
  acf <- .format.acf(acf, N)
  p <- .get.p(dmu, dacf)
  Mu <- .format.mu(mu = mu, dmu = dmu, d2mu = d2mu,
                   N = N, p = p, grad.only = FALSE)
  mu <- Mu$mu
  dmu <- Mu$dmu
  d2mu <- Mu$d2mu
  Dacf <- .format.dacf(dacf = dacf, d2acf = d2acf,
                       N = N, p = p, grad.only = FALSE)
  dacf <- Dacf$dacf
  d2acf <- Dacf$d2acf

  Z <- X - mu
  SigZ <- c(solve(acf, Z))
  SigMu <- solve(acf, dmu)
  SigZ2 <- matrix(NA, N, p)
  for(ii in 1:p) {
    SigZ2[,ii] <- toep.mult(dacf[,ii], SigZ)
  }
  SigZ3 <- solve(acf, SigZ2)
  CP1 <- crossprod(SigMu, SigZ2)
  CP2 <- crossprod(dmu, SigMu)
  CP3 <- colSums(d2mu * SigZ)
  hess <- CP3 - CP1 - t(CP1) - CP2
  for(ii in 1:p) {
    for(jj in ii:p) {
      hess[ii,jj] <- hess[ii,jj] - crossprod(SigZ2[,jj], SigZ3[,ii])
      hess[ii,jj] <- hess[ii,jj] + .5 * crossprod(SigZ, toep.mult(d2acf[,ii,jj], SigZ))
      hess[ii,jj] <- hess[ii,jj] - .5 * (acf$traceT2(d2acf[,ii,jj]) -  acf$traceT4(dacf[,jj], dacf[,ii]))
    }
  }
  hess[lower.tri(hess)] <- t(hess)[lower.tri(hess)]
  hess
}
