#' Gradient of a Stationary Gaussian Log-Likelihood
#'
#' Efficient evaluation of log-likelihood gradient for stationary Gaussian data.
#' @note package "SuperGauss" is required
#' @param X \eqn{n \times 1} matrix, follows N(mean, Variance).
#' @param mean \eqn{n} vector or matrix.
#' @param acf \eqn{n} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param dmean \eqn{n \times p} matrix, where \eqn{p} is the number of parameters, each column is the partial derivative of mean.
#' @param dacf \eqn{n \times p} matrix, each column is the partial derivative of acf.
#' @return The gradient of the log-likelihood.
#' @examples 
#' N <- 30
#' p <- 4
#' X <- as.matrix(rnorm(N))
#' mean <- rnorm(N)
#' acf <- fbm.acf(alpha = 0.8, dT = 1/60, N = N)
#' dmean <- matrix(rnorm(N), N, p)
#' dacf <- matrix(rnorm(N), N, p)
#' acf <- Toeplitz(acf)
#' Snorm.grad(X, mean, acf, dmean, dacf)
#' @export
Snorm.grad <- function(X, mean, acf, dmean, dacf){
  if(is.vector(X)){
    n <- length(X)
    X <- as.matrix(X)
  } else{
    if(ncol(X) != 1){
      stop("X should have only 1 column")
    }
    n <- nrow(X)
  }
  
  if(missing(mean)){
    mean <- rep(0, n)
  } else{
    if(length(mean) == 1){
      mean <- rep(mean, n)
    } else{
      if(length(mean) != n){
        stop("mean has incompatible dimension with X")
      }
      if(is.matrix(mean)){
        mean <- as.vector(mean)
      }
    }
  }
  
  if(class(acf) == "Toeplitz_Matrix"){
    # is Toeplitz
    if(ncol(acf) != n){
      stop("acf has incompatible dimension with X")
    }
  }else{
    if(is.vector(acf)){
      if(length(acf) != n){
        stop("acf has incompatible dimension with X")
      }else{
        acf <- Toeplitz(acf)
      }
    }else{
      stop("acf should be either vector or Toeplitz class")
    }
  }
  
  if(is.vector(dacf)){
    p <- 1
    dacf <- matrix(dacf, ncol = 1)
  } else{
    p <- ncol(dacf)
  }
  if(nrow(dacf) != n){
    stop("dacf has incompatible dimensions with X")
  }
  
  if(missing(dmean)){
    dmean <- matrix(0, n, p)
  } else{
    if(length(dmean) == 1 || length(dmean) == n){
      dmean <- matrix(dmean, n, p)
    } else{
      if(!(is.matrix(dmean) && ncol(dmean) == p && nrow(dmean) == n)){
        stop("dmean has incompatible dimensions with dacf")
      }
    }
  }
  
  X <- X - mean
  SigX <- solve(acf, X)
  trace <- rep(NA, p)
  for(ii in 1:p){
    trace[ii] <- acf$traceT2(dacf[, ii])
  }
  grad <- matrix(NA, p, 1)
  for(ii in 1:p){
    acf$setAcf(dacf[, ii])
    grad[ii] <- crossprod(dmean[, ii], SigX) + crossprod(SigX, acf %*% SigX) / 2
  }
  grad <- grad - trace / 2
  grad
}
