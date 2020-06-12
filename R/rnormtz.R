#' Simulate a stationary Gaussian time series.
#'
#' @param n Number of time series to generate.
#' @param acf Length-`N` vector giving the autocorrelation of the series.
#' @param fft Logical; whether or not to use the `O(N log N)` FFT-based algorithm of Wood and Chan (1994) or the more stable `O(N^2)` Durbin-Levinson algorithm.  See Details.
#' @param Z Optional size `(2N-2) x n` or `N x n` matrix of iid standard normals, to use in the FFT and Durbin-Levinson methods, respectively.
#' @param nkeep Length of time series.  Defaults to `N = length(acf)`. See Details.
#' @param tol Relative tolerance on negative eigenvalues.  See Details.
#'
#' @return Length-`nkeep` vector or size `nkeep x n` matrix with time series as columns.
#'
#' @details The FFT method fails when the embedding circulant matrix is not positive definite.  This is typically due to one of two things:
#' 1.  Roundoff error can make tiny eigenvalues appear negative.  For this purpose, argument `tol` can be used to replace all negative eigenvalues by `tol * ev_max`, where `ev_max` is the largest eigenvalue.
#' 2.  The autocorrelation is decaying too slowly on the given timescale.  To mitigate this, argument `nkeep` can be used to supply a longer `acf` than is required, and keep only the first `nkeep` time series observations.  For consistency, `nkeep` also applies to Durbin-Levinson method.
#'
#' @importFrom fftw FFT IFFT planFFT
#' @example examples/rnormtz.R
#' @export
rnormtz <- function(n = 1, acf, Z, fft = TRUE, nkeep, tol = 1e-6) {
  if(missing(nkeep)) nkeep <- length(acf)
  if(length(acf) <= 2) fft <- FALSE
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
    if(N == 1) {
      X <- sqrt(acf) * Z
    } else {
      X <- DurbinLevinson_ZX(Z = Z, acf = acf)
    }
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
    ## if(N == 0) {
    ##   return(sqrt(acf) * Z[1:(N+1)])
    ## }
    #if(debug) browser()
    # get eigenvalues of circulant embedding
    fft_plan <- planFFT(n = NN)
    psd <- Re(FFT(c(acf, acf[N:2]), plan = fft_plan)[1:(N+1)])
    # clip small negative values
    Smax <- max(psd)
    iS0 <- psd < 0 & abs(psd)/Smax < tol
    psd[iS0] <- tol * Smax
    if(any(psd < 0)) {
      stop("toeplitz(acf) does not embed directly into a positive-definite circulant matrix.\n  Set `fft = FALSE` to use Durbin-Levinson algorithm.")
    }
    psd <- sqrt(psd)
    psd[2:N] <- psd[2:N]/sqrt(2)
    # ifft simulation
    Z <- Z * c(psd, psd[2:N])
    tmp <- Z[2:N,] + 1i * Z[(N+2):NN,]
    Z[NN:(N+2),] <- Z[2:N,] - 1i * Z[(N+2):NN,]
    Z[2:N,] <- tmp
    Z <- apply(Z, 2, IFFT, plan = fft_plan, scale = FALSE)
    # remove roundoff error and discard padded values
    #X <- Re(Z[1:(N+1-ncut),])/sqrt(NN)
    X <- Re(Z[1:nkeep,])/sqrt(NN)
    #if(n != 1) X <- t(X)
  }
  X
}
