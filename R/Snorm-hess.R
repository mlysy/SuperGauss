#' Hessian of a Stationary Gaussian Log-Likelihood
#'
#' Efficient calculation of the Hessian matrix of the log-likelihood of stationary Gaussian data.
#' @note package "SuperGauss" is required
#' @param X \eqn{n \times d} matrix, d i.i.d. vector follows N(mean, Variance).
#' @param mean \eqn{n} vector or matrix.
#' @param acf \eqn{n} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param dmean \eqn{n \times p} matrix, where p is the number of parameters, each column is the partial derivative of mean
#' @param dacf \eqn{n \times p} matrix, each column is the partial deruvative of acf
#' @param d2mean \eqn{n \times p \times p} array
#' @param d2acf \eqn{n \times p \times p} array
#' @note if d2mean and d2acf is \eqn{n \times 1} matrix or array, it only works when p = 1
#' @return The Hessian matrix of the log-likelihood.
#' @examples
#' N <- 30
#' p <- 4
#' X <- as.matrix(rnorm(N))
#' mean <- as.matrix(rnorm(N))
#' acf <- fbm.acf(alpha = 0.8, dT = 1/60, N = N)
#' acf <- Toeplitz(acf = acf)
#' dmean <- matrix(rnorm(N*p), N, p)
#' dacf <- matrix(rnorm(N*p), N, p)
#' d2mean <- array(rnorm(N*p*p), dim = c(N, p, p))
#' d2acf <- array(rnorm(N*p*p), dim = c(N, p, p))
#' Snorm.Hess(X, mean, acf, dmean, dacf, d2mean, d2acf)
#' @export
Snorm.Hess <- function(X, mean, acf, dmean, dacf, d2mean, d2acf){
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
    dacf <- as.matrix(dacf)
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

  if(missing(d2mean)){
    d2mean <- array(0, dim = c(n, p, p))
  } else{
    if(is.vector(d2mean) || is.matrix(d2mean)){
      if(length(d2mean) == 1){
        d2mean <- array(d2mean, dim = c(n, p, p))
      }else{
        if(p == 1 && length(d2mean) == n){
          d2mean <- array(d2mean, dim = c(n, p, p))
        } else{
          stop("d2mean is incompatible with dacf")
        }
      }
    } else{
      if(!prod(as.numeric(dim(d2mean) == c(n, p, p)))){
        stop("dimension of d2mean is incompatible with dacf")
      }
    }
  }

  if(is.vector(d2acf) || is.matrix(d2acf)){
    if(p == 1 && length(d2acf) == n){
      d2acf <- array(d2acf, dim = c(n, p, p))
    } else{
      stop("d2acf is incompatible with dacf")
    }
  } else{
    if(!prod(as.numeric(dim(d2acf) == c(n, p, p)))){
      stop("dimension of d2acf is incompatible with dacf")
    }
  }
  acf.vec <- acf$getAcf()
  X <- X - mean
  SigX <- solve(acf, X)     # stores Sigma^-1 * X
  hess <- matrix(NA, p, p)
  SigMu <- matrix(NA, n, p) # stores Sigma^-1 * Mean_i
  for(ii in 1:p){
    SigMu[,ii] <- solve(acf, dmean[, ii])
  }
  Sig2X <- matrix(NA, n, p) # stores Sigma_i * SigX
  for(ii in 1:p){
    acf$setAcf(dacf[, ii])
    Sig2X[, ii] <- acf %*% SigX
  }
  Sigd2X <- matrix(NA, p, p) # stores SigX' * Sigma_ij * SigX
  for(ii in 1:p){
    for(jj in ii:p){ # symmetric hessian matrix
      acf$setAcf(d2acf[, ii, jj])
      Sigd2X[ii, jj] <- crossprod(SigX, acf %*% SigX)
    }
  }
  acf$setAcf(acf.vec)
  for(ii in 1:p){
    for(jj in ii:p){ # symmetric hessian matrix
      hess[ii, jj] <- crossprod(d2mean[, ii, jj], SigX) - crossprod(SigMu[,ii], Sig2X[, jj]) -
        crossprod(SigMu[,jj], Sig2X[, ii]) - crossprod(dmean[, ii], solve(acf, dmean[, jj])) -
        crossprod(Sig2X[, jj], solve(acf, Sig2X[, ii])) - acf$traceT2(d2acf[, ii, jj]) / 2 +
        acf$traceT4(dacf[, jj], dacf[, ii]) / 2
    }
  }
  hess <- hess + Sigd2X / 2
  hess[lower.tri(hess)] <- t(hess)[lower.tri(hess)]
  hess
}
