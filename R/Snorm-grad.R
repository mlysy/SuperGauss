#' @title Gradient of a stationary Gaussian log-Likelihood
#'
#' @description Efficient evaluation of log-likelihood gradient for stationary Gaussian data.
#' @param X size \code{N x d} matrix, each column i.i.d. follows multivariate Gaussian distribution with mean \code{mu} and Toeplitz variance given by \code{acf}
#' @param mu length \code{N} vector or matrix
#' @param acf length \code{N} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param dmu size \code{N x p} matrix, where \code{p} is the number of parameters, each column is the partial derivative of \code{mu}
#' @param dacf size \code{N x p} matrix, each column is the partial derivative of \code{acf}
#' @note 
#' Order of partial derivative in \code{dmu} and \code{dacf} must be identical. For i-th parameter \eqn{\theta_i}, i-th
#' column of \code{dmu} must be \eqn{\frac{\partial \mu}{\partial \theta_i}}{d\mu/d\theta_i}
#' and i-th column of \code{dacf} must be 
#' \eqn{\frac{\partial acf}{\partial \theta_i}}{dacf/d\theta_i}
#' @return Length \code{N} vector containing the gradient of the log-likelihood.
#' @examples 
#' N <- 300
#' d <- 4
#' X <- matrix(rnorm(N*d), N, d)
#' theta <- 0.1
#' lambda <- 2
#' 
#' mu <- theta^2 * rep(1, N)
#' acf <- exp(-lambda * (1:N - 1))
#' acf <- Toeplitz(acf = acf)
#' dmu <- matrix(0, N, 2)
#' dmu[, 1] <- 2 * theta * rep(1, N)
#' dacf <- matrix(0, N, 2)
#' dacf[, 2] <- -lambda * exp(-lambda * (1:N - 1))
#' 
#' Snorm.grad(X, mu, acf, dmu, dacf)
#' @export
Snorm.grad <- function(X, mu, acf, dmu, dacf){
  if(is.vector(X)){
    N <- length(X)
    d <- 1
    X <- as.matrix(X)
  } else{
    N <- nrow(X)
    d <- ncol(X)
  }
  
  mu <- .format.mu(mu, N)
  acf <- .format.acf(acf, N)
  dacf <- .format.dacf(dacf, N)$dacf
  p <- .format.dacf(dacf, N)$p
  dmu <- .format.dmu(dmu, N, p)
  
  X <- X - mu
  SigX <- solve(acf, X)
  grad <- rep(NA, p)
  for(ii in 1:p){
    grad[ii] <- .trace(crossprod(SigX, toep.mult(dacf[, ii], SigX))) - d * acf$traceT2(dacf[, ii])
  }
  grad <- grad / 2
  grad <- grad + apply(crossprod(dmu, SigX), 1, sum)
  grad
}
