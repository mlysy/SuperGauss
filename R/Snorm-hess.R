#' @title Hessian of a stationary Gaussian log-Likelihood
#'
#' @description Efficient calculation of the Hessian matrix of the log-likelihood of stationary Gaussian data.
#' @param X size \code{N x d} matrix, each column i.i.d. follows multivariate Gaussian distribution with mean \code{mu} and Toeplitz variance given by \code{acf}
#' @param mu length \code{N} vector or matrix
#' @param acf length \code{N} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param dmu size \code{N x p} matrix, where \code{p} is the number of parameters, each column is the partial derivative of \code{mu}
#' @param dacf size \code{N x p} matrix, each column is the partial derivative of \code{acf}
#' @param d2mu size \code{N x p x p} array, each column is the second partial derivative of \code{mu}
#' @param d2acf size \code{N x p x p} array,  each column is the second partial derivative of \code{mu}
#' @note 
#' Order of partial derivative in \code{dmu}, \code{dacf}, \code{d2mu} and \code{d2acf} must be identical. 
#' For i-th and j-th parameter \eqn{\theta_i, \theta_j}, ij-th
#' column of \code{d2mu} must be \eqn{\frac{\partial^2 \mu}{\partial \theta_i \partial \theta_j}}{d2\mu/(d\theta_i d\theta_j)}
#' and ij-th column of \code{d2acf} must be 
#' \eqn{\frac{\partial^2 acf}{\partial \theta_i \partial \theta_j}}{d2acf/(d\theta_i d\theta_j)}
#' @return Size \code{p x p} Hessian matrix of the log-likelihood.
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
#' dmu <- dacf <- matrix(0, N, 2)
#' dmu[, 1] <- 2 * theta * rep(1, N)
#' dacf[, 2] <- -lambda * exp(-lambda * (1:N - 1))
#' d2mu <- d2acf <- array(0, c(N, 2, 2))
#' d2mu[, 1, 1] <- 2 * rep(1, N)
#' d2acf[, 2, 2] <- lambda^2 * exp(-lambda * (1:N - 1))
#' 
#' Snorm.Hess(X, mu, acf, dmu, dacf, d2mu, d2acf)
#' @export
Snorm.Hess <- function(X, mu, acf, dmu, dacf, d2mu, d2acf){
  if(is.vector(X)){
    N <- length(X)
    p <- 1
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
  d2mu <- .format.d2mu(d2mu, N, p)
  d2acf <- .format.d2acf(d2acf, N, p)
  
  X <- X - mu
  SigX <- solve(acf, X)     # stores Sigma^-1 * X, size N x d
  hess <- matrix(NA, p, p)
  SigMu <- solve(acf, dmu) # stores Sigma^-1 * mu_i, size N x p
  Sig2X <- matrix(NA, N, p * d) # stores Sigma_i * SigX, size N x pd
  for(ii in 1:p){
    Sig2X[, d * (ii-1) + 1:d] <- toep.mult(dacf[, ii], SigX)
  }
  for(ii in 1:p){
    for(jj in ii:p){
      hess[ii, jj] <- sum(crossprod(d2mu[, ii, jj], SigX) - crossprod(SigMu[, ii], Sig2X[, d * (jj-1) + 1:d]) - 
                            crossprod(SigMu[, jj], Sig2X[, d * (ii-1) + 1:d])) - d * crossprod(dmu[, ii], SigMu[, jj])
      hess[ii, jj] <- hess[ii, jj] - .trace(crossprod(Sig2X[, d * (jj-1) + 1:d], solve(acf, Sig2X[, d * (ii-1) + 1:d])) - 
                                              crossprod(SigX, toep.mult(d2acf[, ii, jj], SigX)) / 2)
      hess[ii, jj] <- hess[ii, jj] - d / 2 * (acf$traceT2(d2acf[, ii, jj]) - acf$traceT4(dacf[, jj], dacf[, ii]))
    }
  }
  hess[lower.tri(hess)] <- t(hess)[lower.tri(hess)]
  hess
}
