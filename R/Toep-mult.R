#' Toeplitz matrix multiplication.
#'
#' Efficient matrix multiplication with Toeplitz matrix and arbitrary matrix or vector.
#'
#' @param acf Length-\code{N} vector giving the first column (or row) of the Toeplitz matrix.
#' @param X Vector or matrix of compatible dimensions with \code{acf}.
#' @return An \code{N}-row matrix corresponding to \code{toeplitz(acf) \%*\% X}.
#' @importFrom fftw FFT IFFT planFFT
#' @examples
#' N <- 20
#' d <- 3
#' acf <- exp(-(1:N))
#' X <- matrix(rnorm(N*d), N, d)
#' toep.mult(acf, X)
#' @export
toep.mult <- function(acf, X){
  N <- length(acf)
  if(is.vector(X)) X <- matrix(X)
  if(nrow(X) != N) stop("X and acf have incompatible dimensions.")
  d <- ncol(X)

  ## if(is.vector(X)){
  ##   if(length(X) != N){
  ##     stop("X has incompatible dimension with Toeplitz matrix")
  ##   }else{
  ##     d <- 1
  ##   }
  ## }else{
  ##   if(is.matrix(X)){
  ##     if(nrow(X) != N){
  ##       stop("X has incompatible dimension with Toeplitz matrix")
  ##     }else{
  ##       d <- ncol(X)
  ##     }
  ##   }else{
  ##     stop("X should be either vector or matrix")
  ##   }
  ## }

  Y <- matrix(NA, N, d)
  planFFT <- planFFT(n = 2*N-2)
  for(ii in 1:d) {
    vecTemp <- FFT(x = c(acf, rev(acf[-c(1,N)])), plan = planFFT)
    vecTemp <- vecTemp * FFT(x = c(X[,ii], rep(0, N-2)), plan = planFFT)
    Y[,ii] <- IFFT(x = vecTemp, plan = planFFT)[1:N]
  }
  Re(Y)
}
