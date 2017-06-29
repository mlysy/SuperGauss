#' @title Toeplitz matrix multiplication
#' 
#' @description Efficient computation of symmetric Toeplitz matrix times arbitary vector or matrix using Fast Fourier Transformation.
#' @param acf length \code{N} vector giving the first column (or row) of the Toeplitz matrix
#' @param X size \code{N x d} matrix
#' @return Size \code{N x d} matrix.
#' @details Computes \code{toeplitz(acf) * X}.
#' @importFrom fftw FFT IFFT planFFT
#' @examples 
#' N <- 20
#' d <- 3
#' acf <- exp(-(1:N - 1))
#' X <- matrix(rnorm(N*d), N, d)
#' toep.mult(acf, X)
#' @export 
toep.mult <- function(acf, X){
  N <- length(acf)
  
  if(is.vector(X)){
    if(length(X) != N){
      stop("X has incompatible dimension with Toeplitz matrix")
    }else{
      d <- 1
    }
  }else{
    if(is.matrix(X)){
      if(nrow(X) != N){
        stop("X has incompatible dimension with Toeplitz matrix")
      }else{
        d <- ncol(X)
      }
    }else{
      stop("X should be either vector or matrix")
    }
  }
  
  Y <- matrix(NA, N, d)
  planFFT <- planFFT(n = 2 * N)
  for(ii in 1:d){
    vecTemp <- FFT(x = c(acf, 0, rev(acf[-1])), plan = planFFT)
    vecTemp <- vecTemp * FFT(x = c(X[, ii], rep(0, N)), plan = planFFT)
    Y[, ii] <- IFFT(x = vecTemp, plan = planFFT)[1:N]
  }
  Re(Y)
}
