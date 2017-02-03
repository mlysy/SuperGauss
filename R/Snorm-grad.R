#' Gradient of a Stationary Gaussian Log-Likelihood
#'
#' Efficient evaluation of log-likelihood gradient for stationary Gaussian data.
#' @note package "SuperGauss" is required
#' @param X \code{n x d} matrix, d i.i.d. vector follows N(mean, Variance).
#' @param mean \code{n} vector or matrix.
#' @param acf \code{n} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param dmean \code{n x p} matrix, where \code{p} is the number of parameters, each column is the partial derivative of mean.
#' @param dacf \code{n x p} matrix, each column is the partial deruvative of acf.
#' @return The gradient of the log-likelihood.
#' @export
Snorm.grad <- function(X, mean, acf, dmean, dacf, Toep, debug = FALSE){
  if(debug){
    browser()
  }
  if(is.vector(X)){
    n <- length(X)
    X <- matrix(X, n, 1)
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
  
  if(length(acf) != n){
    stop("acf has incompatible dimension with X")
  }
  
  if(missing(Toep)){
    Toep <- new(Toeplitz, n)
  } else{
    if(class(Toep) != "Rcpp_Toeplitz"){
      stop("Toep should be of class \"Toeplitz\"")
    } else{
      if(Toep$DimCheck() != n){
        stop("Toep has incompatible dimension with X")
      }
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
      stop("dmean has incompatible dimensions with dacf")
    }
  }
  
  Toep$AcfInput(acf)
  X <- X - mean
  SigX <- Toep$Solve(X)
  trace <- rep(NA, p)
  for(ii in 1:p){
    trace[ii] <- Toep$TraceProd(dacf[, ii])
  }
  grad <- matrix(NA, p, 1)
  for(ii in 1:p){
    Toep$AcfInput(dacf[, ii])
    grad[ii] <- crossprod(dmean[, ii], SigX) + crossprod(SigX, Toep$Mult(SigX)) / 2
  }
  grad <- grad - trace / 2
  grad
}
