#' @name Choleski
#' @aliases cholXZ cholZX
#' @title Choleski multiplication with Toeplitz variance matrices.
#' @description Multiplies the Choleski decomposition of the Toeplitz matrix with another matrix, or solves a system of equations with the Cholesky factor.
#' @param X Length-\code{N} or \code{N x p} matrix of observations.
#' @param Z Length-\code{N} or \code{N x p} matrix of residuals.
#' @param acf Length-\code{N} autocorrelation vector of the Toeplitz variance matrix.
#' @return Size \code{N x p} residual or observation matrix.
#' @details If \code{C == t(chol(toeplitz(acf)))}, then \code{cholZX} computes \code{C \%*\% Z} and \code{cholZX} computes \code{solve(C, X)}.  Both functions use the Durbin-Levinson algorithm.
#' @rdname Choleski
#' @examples
#' N <- 10
#' p <- 2
#' W <- matrix(rnorm(N * p), N, p)
#' acf <- exp(-(1:N - 1))
#' cholZX(Z = W, acf = acf)
#' cholXZ(X = W, acf = acf)
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

