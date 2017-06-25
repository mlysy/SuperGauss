#' Simulation of a Stationary Gaussian Process.
#'
#' Efficient simulation of a stationary Gaussian process from a given autocorrelation function.
#' @param n number of time series to generate.
#' @param acf vector of length \code{N} giving the autocorrelation of the series.
#' @param fft logical, whether or not to use the superfast fft-based algorithm of Dietrich-Newsam or the more stable Durbin-Levinson algorithm.  See Details.
#' @param Z optional \code{(2*N-2) x n} or \code{N x n} matrix of iid standard normals, for the fft and Durbin-Levinson methods respectively.
#' @param nkeep length of time series.  Defaults to \code{N = length(acf)}.  See Details.
#' @param tol relative tolerance on negative eigenvalues.  See Details.
#' @details The superfast method fails when the circulant matrix is not positive definite.  This is typically due to one of two things, for which the fft algorithm can be tuned with \code{tol} and \code{nkeep}:
#' \itemize{
#' \item \code{tol}: Roundoff error can make tiny eigenvalues appear negative. If \code{evMax} is the maximum eigenvalue, then all negative eigenvalues of magnitude less than \code{tol * evMax} are mapped to this threshold value.  This does not guarantee a positive definite embedding.
#' \item \code{nkeep}: The autocorrelation is decaying too slowly on the given timescale. To mitigate this it is possible to increase the time horizon, i.e. input a longer \code{acf} and keep the first \code{nkeep} time series observations. For consistency, \code{nkeep} also applies to Durbin-Levinson method.
#' }
#' @return A \code{nkeep x n} matrix with time series as columns.
#' @importFrom fftw FFT IFFT planFFT
#' @examples
#' acf <- fbm.acf(alpha = 0.8, dT = 1/60, N = 200)
#' rSnorm(n = 3, acf = acf)
#' @export
rSnorm <- function(n, acf, Z, fft = TRUE, ncut = 0, tol = 1e-6) {
  if(!fft) {
    if(missing(Z)) {
      N <- length(acf)
      Z <- matrix(rnorm(N*n), N, n)
    } else {
      Z <- as.matrix(Z)
      if(length(Z) != n*N) {
        stop("Z has incompatible dimensions with n and acf.")
      }
    }
    X <- DurbinLevinson_ZX(Z = Z, acf = acf)
    if(n == 1) X <- c(X)
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
    X <- Re(Z[1:(N+1-ncut),])/sqrt(NN)
    #if(n != 1) X <- t(X)
  }
  X
}
