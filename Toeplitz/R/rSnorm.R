#' FFT-Based Simulation of Stationary Gaussian Process.
#'
#' @param Z optional \code{(2*N-2) x n} matrix of iid standard normals to use in circulant embedding.
#' @param fft logical, whether or not to use the fft-based algorithm of Dietrich & Newsam (1997) or that of Durbin-Levinson.  The former is much faster but doesn't apply to all Toeplitz variance matrices.
#' @param fft.plan optional plan to pass to fftw.
#' @param ncut length of time series to cut from the end for fft method.  This way, \code{acf} can be "padded" to embed into a circulant matrix and then cut to the correct size.  See details.
#' @param tol relative tolerance on negative eigenvalues.  See details.
#' @details TODO: Explain exactly what ncut and tol do.  Basically, acf can be made longer than needed, and ncut values of the time series are cut off the end.  tol is used to avoid roundoff error on tiny eigenvalues which makes them look negative, i.e., all negative eigenvalues which are less than a \code{tol} fraction of the maximum eigenvalue \code{evMax} are set to \code{tol * evMax}.
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
      stop("toeplitz(acf) does not embed directly into a +ve definite circulant matrix.  Set fft = FALSE to use Durbin-Levinson algorithm.")
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

#' Toeplitz matrix multiplication
#'
#' @details requires packages \code{fftw}.
Toep.mult <- function(acf, x){
  M <- ncol(x)
  N <- nrow(x)
  if(N != length(acf)){
    stop("acf and x have incompatible dimensions.")
  }
  fft.plan <- planFFT(2*N)
  #
  rst <- matrix(NA, N, M)
  acf.fft <- FFT(x = c(acf, 0, acf[N:2]), plan = fft.plan)
  rst <- sapply(1:M, function(ii){
    Re(IFFT(acf.fft * FFT(x = c(x[,ii], rep(0, N)), plan = fft.plan), plan = fft.plan))[1:N]
  })
  rst
}


#' density function of multivariant Normal distribution
#' @note package "Toeplitz" is required
#' @param X \code{M x n} matrix
#' @param mu \code{n} vector
#' @param acf \code{n} vector, or a envirnnment, Toeplitz class
#' @param log return the log-density of True
#' @export
dsNorm <- function(X, mu, acf, log = FALSE){
  M <- ncol(X)
  N <- nrow(X)
  if(length(mu) != N){
    stop("X has incompatible dimensions with mu.")
  }
  if(!(is.vector(acf)||(is.environment(acf)))){
    stop("data type of acf should be vector or Toeplitz class")
  }
  if(is.vector(acf)){
    if(length(acf) != N){
      stop("X has incompatible dimensions with acf.")
    }
    temp <- acf
    acf <- new(Toeplitz, N, M)
    invisible(acf$Compute(temp))
  }
  density <- diag(t(X - mu) %*% acf$InverseProd(X-mu)) + N * log(2*pi) + acf$Det()
  density <- density / -2
  if(log){
    return(density)
  }
  return(exp(density))
}

#' score function of multivariant Normal distribution
#' @note package "Toeplitz" is required
#' @param X \code{n} vector
#' @param mu \code{n} vector
#' @param acf \code{n} vector, or a envirnnment, Toeplitz class
#' @param dmu \code{n x p} matrix, where p is the number of parameters in \code{Theta}
#' @param dacf \code{n x p} matrix
#' @export
sNorm.score <- function(X, mu, acf, dmu, dacf){
  N <- nrow(X)
  p <- ncol(dmu)
  if(length(mu) != N){
    stop("X has incompatible dimensions with mu.")
  }
  if((nrow(dmu) != N)||(nrow(dacf) != N)||(ncol(dacf) != p)){
    stop("dmu and dacf has incompatible dimensions.")
  }
  if(!(is.vector(acf)||(is.environment(acf)))){
    stop("data type of acf should be vector or Toeplitz class")
  }
  if(is.vector(acf)){
    if(length(acf) != N){
      stop("X has incompatible dimensions with acf.")
    }
    temp <- acf
    acf <- new(Toeplitz, N, 1)
    invisible(acf$Compute(temp))
  }
  Xprod <- acf$InverseProd(X-mu)
  score <- sapply(1:p, function(ii){
    -t(dmu[,ii]) %*% Xprod + 1/2 * t(Xprod) %*% Toep.mult(dacf[,ii], Xprod) - 1/2 * acf$TraceProd(dacf[,ii])
  })
  return(score)
} 

#' Hessian matrix of multivariant Normal distribution
#' @note package "Toeplitz" is required
#' @param X \code{n} vector
#' @param mu \code{n} vector
#' @param acf \code{n} vector, or a envirnnment, Toeplitz class
#' @param dmu \code{n x p} matrix, where p is the number of parameters in \code{Theta}
#' @param dacf \code{n x p} matrix
#' @param d2mu \code{n x p x p} array
#' @param d2acf \code{n x p x p} array
#' @export
sNorm.Hess <- function(X, mu, acf, dmu, dacf, d2mu, d2acf){
  N <- length(X)
  p <- ncol(dmu)  
  {
    "dimension check for mu, acf, dmu, dacf, d2mu, d2acf"
    "make sure that d2mu d2acf are array data"
  }
  if(length(mu) != N){
    stop("X has incompatible dimensions with mu.")
  }
  if(length(acf) != N){
    stop("X has incompatible dimensions with acf")
  }
  if(nrow(dmu) != N){
    stop("X has incompatible dimensions with dmu")
  }
  if(nrow(dacf) != N || ncol(dacf) != p){
    stop("dmu has incompatible dimensions with dacf")
  }
  if(!is.array(d2mu) || !is.array(d2acf)){
    stop("d2mu and d2acf should be array data")
  }
  if(!prod(as.numeric(dim(d2mu) == c(N, p, p)))){
    stop("dimension of d2mu is incompatible with X and dmu")
  }
  if(!prod(as.numeric(dim(d2acf) == c(N, p, p)))){
    stop("dimension of d2acf is incompatible with X and dmu")
  }
  temp <- acf
  acf <- new(Toeplitz, N, 1)
  Comp <- (acf$Compute(temp))
  ## in obtaining Hessian_{i,j}
  phi <- Comp$phi *  Comp$sigma2
  dphi[,jj] <- -acf$InverseProd(Toep.mult(dacf[,jj], phi))
  
  #
  Xprod <- acf$InverseProd(X - mu)
  Mprod <- acf$InverseProd(dmu)
  Hess[ii,jj] <- -t(d2mu[ii,jj,]) %*% Xprod + t(Mprod[,ii]) %*% Toep.mult(dacf[,jj],Xprod) -
                 t(Mprod[,jj]) %*% Toep.mult(dacf[,ii],Xprod) + t(dmu[,ii]) %*% Mprod[,jj] +
                 t(Xprod) %*% Toep.mult(dacf[,jj], acf$InverseProd(Toep.mult(dacf[, ii], Xprod))) + 
                 1/2 * t(Xprod) %*% Toep.mult(d2acf[ii,jj,], Xprod) - 1/2 * "trace part"
}
