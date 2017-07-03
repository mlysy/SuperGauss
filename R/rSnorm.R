#' Simulation of a stationary Gaussian time series.
#'
#' @param n Number of time series to generate.
#' @param acf Length-\code{N} vector giving the autocorrelation of the series.
#' @param fft Logical, whether or not to use the superfast FFT-based algorithm of Chan and Wood or the more stable Durbin-Levinson algorithm.  See Details.
#' @param Z Optional size \code{(2N-2) x n} or \code{N x n} matrix of iid standard normals, to use in the fft and Durbin-Levinson methods respectively.
#' @param nkeep Length of time series.  Defaults to \code{N = length(acf)}. See Details.
#' @param tol Relative tolerance on negative eigenvalues.  See Details.
#' @details The superfast method fails when the circulant matrix is not positive definite.  This is typically due to one of two things, for which the FFT algorithm can be tuned with \code{tol} and \code{nkeep}:
#' \describe{
#' \item{\code{tol}}{Roundoff error can make tiny eigenvalues appear negative. If \code{evMax} is the maximum eigenvalue, then all negative eigenvalues of magnitude less than \code{tol * evMax} are mapped to this threshold value.  This does not guarantee a positive definite embedding.}
#' \item{\code{nkeep}}{The autocorrelation is decaying too slowly on the given timescale. To mitigate this it is possible to increase the time horizon, i.e. input a longer \code{acf} and keep the first \code{nkeep} time series observations. For consistency, \code{nkeep} also applies to Durbin-Levinson method.}
#' }
#' @return Length-\code{nkeep} vector or size \code{nkeep x n} matrix with time series as columns.
#' @importFrom fftw FFT IFFT planFFT
#' @examples
#' N <- 10
#' acf <- exp(-(1:N - 1))
#' rSnorm(n = 3, acf = acf)
#' @export
rSnorm <- function(n = 1, acf, Z, fft = TRUE, nkeep, tol = 1e-6) {
  if(missing(nkeep)) nkeep <- length(acf)
  if(!fft) {
    if(missing(Z)) {
      N <- length(acf)
      Z <- matrix(rnorm(N*n), N, n)
    } else {
      Z <- as.matrix(Z)
      if(length(Z) != n*N) {
        stop("Z has incompatible dimensions with n and length(acf).")
      }
    }
    X <- DurbinLevinson_ZX(Z = Z, acf = acf)
    X <- X[1:nkeep,]
    #if(n == 1) X <- c(X)
  } else {
    N <- length(acf)-1
    NN <- max(2*N,1)
    if(missing(Z)) {
      Z <- matrix(rnorm(n*NN),NN,n)
    } else {
      Z <- as.matrix(Z)
      if(length(Z) != n*NN) {
        stop("Z has incompatible dimensions with n and acf.")
      }
    }
    if(N == 0) {
      return(sqrt(acf) * Z[1:(N+1)])
    }
    #if(debug) browser()
    # get eigenvalues of circulant embedding
    fft.plan <- planFFT(n = NN)
    psd <- Re(FFT(c(acf, acf[N:2]), plan = fft.plan)[1:(N+1)])
    # clip small negative values
    Smax <- max(psd)
    iS0 <- psd < 0 & abs(psd)/Smax < tol
    psd[iS0] <- tol * Smax
    if(any(psd < 0)) {
      stop("toeplitz(acf) does not embed directly into a +ve definite circulant matrix.\n  Set fft = FALSE to use Durbin-Levinson algorithm.")
    }
    psd <- sqrt(psd)
    psd[2:N] <- psd[2:N]/sqrt(2)
    # ifft simulation
    Z <- Z * c(psd, psd[2:N])
    tmp <- Z[2:N,] + 1i * Z[(N+2):NN,]
    Z[NN:(N+2),] <- Z[2:N,] - 1i * Z[(N+2):NN,]
    Z[2:N,] <- tmp
    Z <- apply(Z, 2, IFFT, plan = fft.plan, scale = FALSE)
    # remove roundoff error and discard padded values
    #X <- Re(Z[1:(N+1-ncut),])/sqrt(NN)
    X <- Re(Z[1:nkeep,])/sqrt(NN)
    #if(n != 1) X <- t(X)
  }
  X
}
