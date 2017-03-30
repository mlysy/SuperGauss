#' @name Choleski
#' @aliases cholXZ cholZX
#' @title Toeplitz variance matrix Choleski decomposition
#' @description Compute Choleski decomposition and multiply by another matrix, or solve Choleski system of equations.
#' @param X \code{n x p} matrix of observations.
#' @param Z \code{n x p} matrix of residuals.
#' @param acf vector of length \code{n}, autocorrelation of Toeplitz matrix.
#' @details \code{cholZX} computes \code{X = chol(toeplitz(acf))' Z}, and \code{cholXZ} computes \code{Z = solve(chol(toeplitz(acf))', X)}.  Both are done with Durbin-Levinson algorithm.
#' @return An \code{n x p} matrix.
#' @rdname Choleski
#' @export
cholZX <- function(Z, acf) {
  n <- length(acf)
  Z <- as.matrix(Z)
  if(nrow(Z) != n) stop("Z and acf have incompatible dimensions.")
  toeplitzZX(Z = Z, acf = acf)
}

#' @rdname Choleski
#' @export
cholXZ <- function(X, acf) {
  n <- length(acf)
  X <- as.matrix(X)
  if(nrow(X) != n) stop("X and acf have incompatible dimensions.")
  toeplitzXZ(X = X, acf = acf)
}

