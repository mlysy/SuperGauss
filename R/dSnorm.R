#' Density of a Stationary Gaussian Process
#'
#' Efficient density evaluation for the multivariate normal distribution with Toeplitz variance matrix.
#' @note package "SuperGauss" is required
#' @note When input {X, mu, acf} data type can be either vector or matrix, if vector, convert it into \code{n x 1} matrix
#' @param X \eqn{N \times d} matrix, d i.i.d. vector follows N(mu, Variance)
#' @param mu \eqn{N} vector or matrix
#' @param acf \eqn{N} vector or matrix, first column of variance matrix, or a Toeplitz class initialized by acf
#' @param log logical, if \code{TRUE} returns log-densities.
#' @param return vector of densities.
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
#' 
#' dSnorm(X, mu = mean, acf = acf, log = TRUE)
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
