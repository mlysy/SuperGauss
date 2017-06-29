#' @title Joint density of a stationary Gaussian process
#'
#' @description  Efficient density evaluation for the multivariate normal distribution with Toeplitz variance matrix.
#' @param X size \code{N x d} matrix, each column i.i.d. follows multivariate Gaussian distribution with mean \code{mu} and Toeplitz variance given by \code{acf}
#' @param mu length \eqn{N} vector or matrix, mean of the process
#' @param acf length \code{N} vector containing the first column of process variance, or a size \code{N} Toeplitz class initialized by acf
#' @param log logical, if \code{TRUE} returns log-densities
#' @return Density (log density) of stationary Gaussian system.
#' @examples 
#' N <- 10
#' d <- 4
#' X <- matrix(rnorm(N*d), N, d)
#' theta <- 0.1
#' lambda <- 2
#' 
#' mu <- theta^2 * rep(1, N)
#' acf <- exp(-lambda * (1:N - 1))
#' acf <- Toeplitz(acf = acf)
#' 
#' dSnorm(X, mu, acf, log = TRUE)
#' @export
dSnorm <- function(X, mu, acf, log = TRUE){
  if(is.vector(X)){
    N <- length(X)
    d <- 1
    X <- as.matrix(X)
  } else{
    N <- nrow(X)
    d <- ncol(X)
  }
  
  # formating mu and acf
  mu <- .format.mu(mu, N)
  acf <- .format.acf(acf, N)
  
  X <- X - mu
  density <- .trace(crossprod(X, solve(acf, X)))
  density <- density + determinant(acf) * d + N * d * log(2 * pi)
  density <- density / -2
  
  if(log){
    density
  }
  else{
    exp(density)
  }
}
