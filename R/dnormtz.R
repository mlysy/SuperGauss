#' Density of a multivariate normal with Toeplitz variance matrix.
#'
#' @param X Vector of length `N` or `N x n` matrix, of which each column is a multivariate observation.
#' @param mu Vector or matrix of mean values of compatible dimensions with `X`.  Defaults to all zeros.
#' @param acf Vector of length `N` containing the first column of the Toeplitz variance matrix.
#' @param log Logical; whether to return the multivariate normal density on the log scale.
#' @param method Which calculation method to use.  Choices are: `gschur` for a modified version of the Generalized Schur algorithm of Ammar & Gragg (1988), or `ltz` for the Levinson-Trench-Zohar method.  The former scales as `O(N log^2 N)` whereas the latter scales as `O(N^2)` and should only be used for `N < 300`.
#' @return Vector of `n` (log-)densities, one for each column of `X`.
#' @example examples/dnormtz.R
#' @export
dnormtz <- function(X, mu, acf, log = FALSE, method = c("gschur", "ltz")) {
  # format arguments
  method <- match.arg(method)
  if(is.vector(X)) X <- as.matrix(X)
  N <- nrow(X)
  d <- ncol(X)
  if(missing(mu)) mu <- matrix(0, N, d)
  Z <- X - mu
  check_acf(acf, N)
  ## acf <- .format.acf(acf, N)
  # log-density calculation
  if(method == "gschur") {
    Tz <- Toeplitz$new(acf = acf)
    IP <- colSums(Z * solve(Tz, Z))
    ldV <- determinant(Tz)
  } else if(method == "ltz") {
    # fixme: use LTZ instead of DL!
    DL <- DurbinLevinson_crossprod(X = Z, Y = Z, acf = acf, calc_mode = 2L)
    IP <- DL$IP
    ldV <- DL$ldV
  }
  ld <- -.5 * (IP + ldV + N * log(2 * pi))
  if(!log){
    ld <- exp(ld)
  }
  ld
}

## #' @rdname dSnorm
## #' @details \code{dSnorm} and \code{dSnormDL} have identical outputs, with the former using the generalized Schur algorithm and the latter, the Durbin-Levinson algorithm, which is more common but slower.  \code{dSnormDL} is provided mainly for speed comparisons.
## #' @export
## dSnormDL <- function(X, mu, acf, log = FALSE) {
##   # format arguments
##   if(is.vector(X)) X <- as.matrix(X)
##   N <- nrow(X)
##   d <- ncol(X)
##   if(missing(mu)) mu <- matrix(0, N, d)
##   Z <- X - mu
##   if(length(acf) != N) {
##     stop("X and acf have incompatible dimensions.")
##   }
##   # density calculation
##   DL <- DurbinLevinson_Eigen(X = Z, Y = Z, acf = acf, calcMode = 2L)
##   ld <- -.5 * (as.numeric(DL$IP) + DL$ldV + N * log(2 * pi))
##   if(!log){
##     ld <- exp(ld)
##   }
##   ld
## }
