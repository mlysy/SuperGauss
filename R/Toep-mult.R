#' Toeplitz matrix multiplies vector or matrix
#' 
#' efficient computation of symmetric Toeplitz matrix times arbitary vector or matrix using Fast Fourier Transformation
#' @param vec vector of length \code{N} giving the first column (or row) of the Toeplitz matrix
#' @param X vector of length \code{N} or matrix of dimension \code{N x d}
#' @return A matrix of size of X
#' @details computes \code{toeplitz(vec) %*% X}
#' computes 
#' @importFrom fftw FFT IFFT planFFT
#' @examples 
#' vec <- rnorm(200)
#' X <- matrix(rnorm(200 * 5), 200, 5)
#' vec.mult(vec, X)
#' @export 
vec.mult <- function(vec, X){
  N <- length(vec)
  
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
    vecTemp <- FFT(x = c(vec, 0, rev(vec[-1])), plan = planFFT)
    vecTemp <- vecTemp * FFT(x = c(X[, ii], rep(0, N)), plan = planFFT)
    Y[, ii] <- IFFT(x = vecTemp, plan = planFFT)[1:N]
  }
  Re(Y)
}
