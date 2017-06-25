#' @name Choleski
#' @aliases cholXZ cholZX
#' @title Toeplitz variance matrix Choleski decomposition
#' @description
#' Compute Choleski decomposition and multiply by another matrix, or solve Choleski system of equations.
#' @param X \eqn{n \times p} matrix of observations.
#' @param Z \eqn{n \times p} matrix of residuals.
#' @param acf vector of length \eqn{n}, autocorrelation of Toeplitz matrix.
#' @details
#' \itemize{
#'   \item{1}{\code{cholZX} computes \code{X = chol(toeplitz(acf))' Z}}
#'   \item{2}{\code{cholXZ} computes\code{Z = solve(chol(toeplitz(acf))', X)}}
#' }
#' Both are done with Durbin-Levinson algorithm.
#' @return An \eqn{n \times p} matrix.
#' @rdname Choleski
#' @examples
#' N <- 30
#' p <- 4
#' Mat <- matrix(rnorm(N * p), N, p)
#' acf <- fbm.acf(alpha = 0.8, dT = 1/60, N = N)
#' cholZX(Z = Mat, acf = acf)
#' cholXZ(X = Mat, acf = acf)
#' @export
cholZX <- function(Z, acf) {
  n <- length(acf)
  Z <- as.matrix(Z)
  if(nrow(Z) != n) stop("Z and acf have incompatible dimensions.")
  DurbinLevinson_ZX(Z = Z, acf = acf)
}

#' @rdname Choleski
#' @export
cholXZ <- function(X, acf) {
  n <- length(acf)
  X <- as.matrix(X)
  if(nrow(X) != n) stop("X and acf have incompatible dimensions.")
  DurbinLevinson_XZ(X = X, acf = acf)
}

