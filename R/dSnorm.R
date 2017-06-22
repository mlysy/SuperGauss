#' Density of a Stationary Gaussian Process
#'
#' Efficient density evaluation for the multivariate normal distribution with Toeplitz variance matrix.
#' @note package "SuperGauss" is required
#' @note When input {X, mean, acf} data type can be either vector or matrix, if vector, convert it into \code{n x 1} matrix
#' @param X \eqn{N \times d} matrix, d i.i.d. vector follows N(mean, Variance)
#' @param mean \eqn{N} vector or matrix
#' @param acf \eqn{N} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param log logical, if \code{TRUE} returns log-densities.
#' @param return vector of densities.
#' @examples 
#' N <- 30
#' d <- 4
#' X <- matrix(rnorm(N*d), N, d)
#' mean <- rnorm(N)
#' acf <- fbm.acf(alpha = 0.8, dT = 1/60, N = N)
#' acf <- Toeplitz(acf)
#' dSnorm(X, mean = mean, acf = acf, log = TRUE)
#' @export
dSnorm <- function(X, mean, acf, log = FALSE){
  if(is.vector(X)){
    n <- length(X)
    X <- as.matrix(X)
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
  
  X <- X - mean
  density <- rep(NA, d)
  for(ii in 1:d){
    density[ii] <- crossprod(X[, ii], solve(acf, X[, ii]))
  }
  density <- density + n * log(2*pi) + determinant(acf)
  density <- density / -2
  if(log){
    density
  }
  else{
    exp(density)
  }
}
