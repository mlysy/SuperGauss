#' Density of a Stationary Gaussian Process
#'
#' Efficient density evaluation for the multivariate normal distribution with Toeplitz variance matrix.
#' @note package "SuperGauss" is required
#' @note When input {X, mean, acf} data type can be either vector or matrix, if vector, convert it into \code{n x 1} matrix
#' @param X \code{n x d} matrix, d i.i.d. vector follows N(mean, Variance)
#' @param mean \code{n} vector or matrix
#' @param acf \code{n} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param log logical, if \code{TRUE} returns log-densities.
#' @param return vector of densities.
#' @export
dSnorm <- function(X, mean, acf, Toep, log = FALSE){
  if(is.vector(X)){
    n <- length(X)
    X <- matrix(X, n, 1)
    d <- 1
  } else{
    n <- nrow(X)
    d <- ncol(X)
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
  
  Toep$AcfInput(acf)
  X <- X - mean
  density <- rep(NA, d)
  for(ii in 1:d){
    density[ii] <- crossprod(X[, ii], Toep$SolveVec(X[, ii]))
  }
  density <- density + n * log(2*pi) + Toep$Det()
  density <- density / -2
  if(log){
    density
  }
  else{
    exp(density)
  }
}
