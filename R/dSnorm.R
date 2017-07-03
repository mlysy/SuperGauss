#' @title Density of a multivariate normal with Toeplitz variance matrix.
#'
#' @description  Efficient density evaluation for the multivariate normal distribution with Toeplitz variance matrix.
#' @param X Vector or matrix, of which each column is a multivariate observation.
#' @param mu Vector or matrix of mean values of compatible dimensions with \code{X}.  Defaults to all zeros.
#' @param acf Vector containing the first column of the Toeplitz variance matrix.  For \code{dSnorm}, can also be a \code{Toeplitz} object.
#' @param log Logical, whether to the return PDF on log scale.
#' @return Vector of (log-)densities, one for each column of \code{X}.
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
dSnorm <- function(X, mu, acf, log = FALSE) {
  # format arguments
  if(is.vector(X)) X <- as.matrix(X)
  N <- nrow(X)
  d <- ncol(X)
  if(missing(mu)) mu <- matrix(0, N, d)
  Z <- X - mu
  acf <- .format.acf(acf, N)
  # log-density calculation
  IP <- colSums(Z * solve(acf, Z))
  ldV <- determinant(acf)
  ld <- -.5 * (IP + ldV + N * log(2 * pi))
  if(!log){
    ld <- exp(ld)
  }
  ld
}

#' @rdname dSnorm
#' @details \code{dSnorm} and \code{dSnormDL} have identical outputs, with the former using the generalized Schur algorithm and the latter, the Durbin-Levinson algorithm, which is more common but slower.  \code{dSnormDL} is provided mainly for speed comparisons.
#' @export
dSnormDL <- function(X, mu, acf, log = FALSE) {
  # format arguments
  if(is.vector(X)) X <- as.matrix(X)
  N <- nrow(X)
  d <- ncol(X)
  if(missing(mu)) mu <- matrix(0, N, d)
  Z <- X - mu
  if(length(acf) != N) {
    stop("X and acf have incompatible dimensions.")
  }
  # density calculation
  DL <- DurbinLevinson_Eigen(X = Z, Y = Z, acf = acf, calcMode = 2L)
  ld <- -.5 * (as.numeric(DL$IP) + DL$ldV + N * log(2 * pi))
  if(!log){
    ld <- exp(ld)
  }
  ld
}
