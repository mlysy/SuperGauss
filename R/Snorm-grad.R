#' Gradient of a Stationary Gaussian Log-Likelihood
#'
#' Efficient evaluation of log-likelihood gradient for stationary Gaussian data.
#' @note package "SuperGauss" is required
#' @param X \eqn{N \times d} matrix, each column i.i.d. follows multivariate Gaussian distribution with mean \code{mu} and Toeplitz variance given by \code{acf}.
#' @param mu length \eqn{N} vector or matrix.
#' @param acf length \eqn{N} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param dmu size \eqn{N \times p} matrix, where \eqn{p} is the number of parameters, each column is the partial derivative of \code{mu}.
#' @param dacf size \eqn{N \times p} matrix, each column is the partial derivative of \code{acf}.
#' @note 
#' the order of partial derivative in \code{dmu} and \code{\dacf} must be identical. Assuming that \eqn{\beta} is
#' the second parameter, then second column of \code{dmu} should be \eqn{\frac{\partial \mu}{\partial \beta}}, second column of 
#' \code{dacf} should be \eqn{\frac{\partial acf}{\partial \beta}}
#' @return The gradient of the log-likelihood.
#' @examples 
#' N <- 300
#' d <- 4
#' X <- matrix(rnorm(N*d), N, d)
#' theta <- 0.1
#' lambda <- 2
#' 
#' mu <- theta^2 * rep(1, N)
#' acf <- exp(-lambda * (1:N - 1))
#' acf <- Toeplitz(acf)
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
    grad[ii] <- .trace(crossprod(SigX, Toep.mult(dacf[, ii], SigX))) - d * acf$traceT2(dacf[, ii])
  }
  grad <- grad / 2
  grad <- grad + apply(crossprod(dmu, SigX), 1, sum)
  grad
}
