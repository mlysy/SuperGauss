#' Hessian of a Stationary Gaussian Log-Likelihood
#'
#' Efficient calculation of the Hessian matrix of the log-likelihood of stationary Gaussian data.
#' @note package "Toeplitz" is required
#' @param X \code{n x d} matrix, d i.i.d. vector follows N(mean, Variance).
#' @param mean \code{n} vector or matrix.
#' @param acf \code{n} vector or matrix, first column of variance matrix
#' @param Toep \code{n x 1} Toeplitz class, space for Toeplitz-related computation.
#' @param dmean \code{n x p} matrix, where p is the number of parameters, each column is the partial derivative of mean
#' @param dacf \code{n x p} matrix, each column is the partial deruvative of acf
#' @param d2mean \code{n x p x p} array
#' @param d2acf \code{n x p x p} array
#' @return The Hessian matrix of the log-likelihood.
#' @export
Snorm.Hess <- function(X, mean, acf, dmean, dacf, d2mean, d2acf, Toep){
  if(is.vector(X)){
    n <- length(X)
    X <- matrix(X, n, 1)
  } else{
    if(ncol(X) != 1){
      stop("X should have only 1 column")
    }
    n <- nrow(X)
    p <- ncol(dmean)
  }
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
  if(length(acf) != n){
    stop("acf has incompatible dimension with X")
  }
  if(missing(Toep)){
    Toep <- new(Toeplitz, n)
  } else{
    if(Toep$DimCheck() != n){
      stop("Toep has incompatible dimension with X")
    }
  }
  if((nrow(dmean) != n)||(nrow(dacf) != n)){
    stop("dmean or dacf has incompatible dimensions with X.")
  }
  if(ncol(dacf) != p){
    stop("dmean has incompatible dimensions with dacf.")
  }
  if(!is.array(d2mean) || !is.array(d2acf)){
    stop("d2mu and d2acf should be array data")
  }
  if(!prod(as.numeric(dim(d2mean) == c(n, p, p)))){
    stop("dimension of d2mu is incompatible with dmean")
  }
  if(!prod(as.numeric(dim(d2acf) == c(n, p, p)))){
    stop("dimension of d2acf is incompatible with dmean")
  }

  X <- X - mean
  Toep$AcfInput(acf)
  SigX <- Toep$Solve(X)     # stores Sigma^-1 * X
  hess <- matrix(NA, p, p)
  SigMu <- matrix(NA, n, p) # stores Sigma^-1 * Mean_i
  for(ii in 1:p){
    SigMu[,ii] <- Toep$SolveVec(dmean[, ii])
  }
  Sig2X <- matrix(NA, n, p) # stores Sigma_i * SigX
  for(ii in 1:p){
    Toep$AcfInput(dacf[, ii])
    Sig2X[, ii] <- Toep$Mult(SigX)
  }
  Sigd2X <- matrix(NA, p, p) # stores SigX' * Sigma_ij * SigX
  for(ii in 1:p){
    for(jj in ii:p){ # symmetric hessian matrix
      Toep$AcfInput(d2acf[, ii, jj])
      Sigd2X[ii, jj] <- crossprod(SigX, Toep$Mult(SigX))
    }
  }
  Toep$AcfInput(acf)
  for(ii in 1:p){
    for(jj in ii:p){ # symmetric hessian matrix
      hess[ii, jj] <- -crossprod(d2mean[, ii, jj], SigX) + crossprod(SigMu[,ii], Sig2X[, jj]) -
        crossprod(SigMu[,jj], Sig2X[, ii]) + crossprod(dmean[, ii], Toep$SolveVec(dmean[, jj])) +
        crossprod(Sig2X[, jj], Toep$SolveVec(Sig2X[, ii])) - Toep$TraceProd(d2acf[, ii, jj]) / 2 +
        Toep$TraceDeriv(dacf[, jj], dacf[, ii]) / 2
    }
  }
  hess <- hess + Sigd2X / 2
  hess[lower.tri(hess)] <- t(hess)[lower.tri(hess)]
  hess
}
