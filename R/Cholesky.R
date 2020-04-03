#' Cholesky multiplication with Toeplitz variance matrices.
#'
#' Multiplies the Cholesky decomposition of the Toeplitz matrix with another matrix, or solves a system of equations with the Cholesky factor.
#'
#' @name Cholesky
#' @aliases cholXZ cholZX
#'
#' @param X Length-`N` or `N x p` matrix of observations.
#' @param Z Length-`N` or `N x p` matrix of residuals.
#' @param acf Length-`N` autocorrelation vector of the Toeplitz variance matrix.
#'
#' @return Size `N x p` residual or observation matrix.
#'
#' @details If `C == t(chol(toeplitz(acf)))`, then `cholZX()` computes `C %*% Z` and `cholZX()` computes `solve(C, X)`.  Both functions use the Durbin-Levinson algorithm.
#' @example examples/Cholesky.R
NULL

#' @rdname Cholesky
#' @export
cholZX <- function(Z, acf) {
  n <- length(acf)
  Z <- as.matrix(Z)
  if(nrow(Z) != n) stop("Z and acf have incompatible dimensions.")
  DurbinLevinson_ZX(Z = Z, acf = acf)
}

#' @rdname Cholesky
#' @export
cholXZ <- function(X, acf) {
  n <- length(acf)
  X <- as.matrix(X)
  if(nrow(X) != n) stop("X and acf have incompatible dimensions.")
  DurbinLevinson_XZ(X = X, acf = acf)
}

