#' Simulation of a Stationary Gaussian Process.
#'
#' Efficient simulation of a stationary Gaussian process from a given autocorrelation function.
#' @param Z optional \code{(2*N-2) x n} matrix of iid standard normals to use in circulant embedding.
#' @param fft logical, whether or not to use a superfast fft-based algorithm or the more stable Durbin-Levinson algorithm.  See details.
#' @param fft.plan optional plan to pass to \code{fftw}.
#' @param ncut length of time series to cut from the end for fft method.  This way, \code{acf} can be "padded" to embed into a circulant matrix and then cut to the correct size.  See details.
#' @param tol relative tolerance on negative eigenvalues.  See details.
#' @details TODO: Explain exactly what ncut and tol do.  Basically, acf can be made longer than needed, and ncut values of the time series are cut off the end.  tol is used to avoid roundoff error on tiny eigenvalues which makes them look negative, i.e., all negative eigenvalues which are less than a \code{tol} fraction of the maximum eigenvalue \code{evMax} are set to \code{tol * evMax}.
#' @importFrom fftw FFT IFFT planFFT
#' @export
rSnorm <- function(n, acf, Z, fft = TRUE, fft.plan,
                   ncut = 0, tol = 1e-6,
                   debug = FALSE) {
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
    X <- toeplitzZX(Z = Z, acf = acf)
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
    if(debug) browser()
    # get eigenvalues of circulant embedding
    if(missing(fft.plan)) fft.plan <- planFFT(n = NN)
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
    if(n != 1) X <- t(X)
  }
  X
}
