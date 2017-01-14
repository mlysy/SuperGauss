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

#' density function of multivariant Normal distribution with specific Toeplitz variance
#' @note package "Toeplitz" is required
#' @note When input {X, mean, acf} data type can be either vector or matrix, if vector, convert it into \code{n x 1} matrix
#' @param X, \code{n x d} matrix, d i.i.d. vector follows N(mean, Variance)
#' @param mean, \code{n} vector or matrix
#' @param acf, \code{n} vector or matrix, first column of variance matrix
#' @param Toep, \code{n x d} Toeplitz class, space for Toeplitz-related computation
#' @param log, logic, return the log-density of True
#' @export
dSnorm <- function(X, mean, acf, Toep, log = FALSE){
  if(is.vector(X)){
    n <- length(X)
    X <- matrix(X, n, 1)
  } else{
    n <- ncol(X)
    d <- nrow(X)
  }
  if(length(mean) == 1){
    mean <- matrix(mean, n, 1)
  } else{
    if(length(mean) != n){
      stop("mean has incompatible dimension with X")
    }
    if(is.vector(mean)){
      mean <- matrix(mean, n, 1)
    }
  }
  if(length(acf) != n){
    stop("acf has incompatible dimension with X")
  }
  if(is.vector(acf)){
    acf <- matrix(acf, n, 1)
  }
  if(missing(Toep)){
    Toep <- new(Toeplitz, n)
  } else{
    if(Toep$DimCheck() != n){
      stop("Toep has incompatible dimension with X")
    }
  }
  Toep$AcfInput(acf)
  X <- X - mean
  density <- crossprod(X, acf$Solve(X))
  density <- density + N * log(2*pi) + Toep$Det()
  density <- density / -2
  if(log){
    density
  }
  else{
    exp(density)
  }
}

#' gradiant function of multivariant Normal distribution
#' @note package "Toeplitz" is required
#' @param X, \code{n x d} matrix, d i.i.d. vector follows N(mean, Variance)
#' @param mean, \code{n} vector or matrix
#' @param acf, \code{n} vector or matrix, first column of variance matrix
#' @param Toep, \code{n x 1} Toeplitz class, space for Toeplitz-related computation
#' @param dmean \code{n x p} matrix, where p is the number of parameters, each column is the partial derivative of mean
#' @param dacf \code{n x p} matrix, each column is the partial deruvative of acf
#' @export
Snorm.grad <- function(X, mean, acf, dmean, dacf, Toep){
  n <- nrow(X)
  p <- ncol(dmean)
  if(length(mean) == 1){
    mean <- rep(mean, n)
  } else{
    if(length(mean) != n){
      stop("mean has incompatible dimension with X")
    } 
  }
  if(length(acf) != n){
    stop("acf has incompatible dimensions with X")
  }
  if((nrow(dmean) != n)||(nrow(dacf) != n)||(ncol(dacf) != p)){
    stop("dmean and dacf has incompatible dimensions.")
  }
  if(missing(Toep)){
    Toep <- new(Toeplitz, n)
  }
  else{
    if(Toep$DimCheck() != n){
      stop("Toep has incompatible dimensions with X")
    }
  }
  X <- X - mean
  Toep$AcfInput(acf)
  SigX <- Toep$Solve(X)
  trace <- rep(NA, p)
  for(ii in 1:p){
    trace[ii] <- Toep$TraceProd(dacf[, ii])
  }
  grad <- rep(NA, p)
  for(ii in 1:p){
    grad.val <- -crossprod(dmean[, ii], SigX)
    Toep$AcfInput(dacf[, ii])
    grad.val <- grad[ii] + crossprod(SigX, Toep$Mult(SigX)) / 2
    grad[ii] <- grad.val
  }
  grad <- grad - trace / 2
  grad
} 

#' Hessian matrix of multivariant Normal distribution
#' @note package "Toeplitz" is required
#' @param X, \code{n x d} matrix, d i.i.d. vector follows N(mean, Variance)
#' @param mean, \code{n} vector or matrix
#' @param acf, \code{n} vector or matrix, first column of variance matrix
#' @param Toep, \code{n x 1} Toeplitz class, space for Toeplitz-related computation
#' @param dmean \code{n x p} matrix, where p is the number of parameters, each column is the partial derivative of mean
#' @param dacf \code{n x p} matrix, each column is the partial deruvative of acf
#' @param d2mean \code{n x p x p} array
#' @param d2acf \code{n x p x p} array
#' @export
Snorm.Hess <- function(X, mean, acf, dmean, dacf, d2mean, d2acf, Toep){
  n <- length(X)
  p <- ncol(dmu)  
  if(length(mean) == 1){
    mean <- rep(mean, n)
  } else{
    if(length(mean) != n){
      stop("mean has incompatible dimension with X")
    } 
  }
  if(length(acf) != n){
    stop("X has incompatible dimensions with acf")
  }
  if(nrow(dmean) != n){
    stop("X has incompatible dimensions with dmean")
  }
  if(nrow(dacf) != n || ncol(dacf) != p){
    stop("dacf has incompatible dimensions with dmean")
  }
  if(!is.array(d2mean) || !is.array(d2acf)){
    stop("d2mu and d2acf should be array data")
  }
  if(!prod(as.numeric(dim(d2mean) == c(n, p, p)))){
    stop("dimension of d2mu is incompatible with X and dmean")
  }
  if(!prod(as.numeric(dim(d2acf) == c(n, p, p)))){
    stop("dimension of d2acf is incompatible with X and dmean")
  }
  if(missing(Toep)){
    Toep <- new(Toeplitz, n)
  } else{
    if(Toep$DimCheck() != n){
      stop("dimension of Toep is incompatible with X")
    }
  }
  X <- X - mean
  Toep$AcfInput(acf)
  SigX <- Toep$Solve(X)     # stores Sigma^-1 * X
  hess <- matrix(NA, p, p)
  SigMu <- matrix(NA, n, p) # stores Sigma^-1 * Mean_i
  for(ii in 1:p){
    SigMu[,ii] <- Toep$Solve(dmean[, ii])
  }
  Sig2X <- matrix(NA, n, p) # stores Sigma_i * SigX
  for(ii in 1:p){
    Toep$AcfInput(dacf[, ii])
    Sig2X <- Toep$Mult(SigX)
  }
  Sigd2X <- matrix(NA, p, p) # stores SigX' * Sigma_ij * SigX
  for(ii in 1:p){
    for(jj in ii:p){ # symmetric hessian matrix
      Toep$AcfInput(d2acf[, ii, jj])
      Sigd2X[ii, jj] <- crossprod(SigX, Toep$Mult(SigX))
    }
  }
  Toep$AcfInput(acf)
  for(ii in 1:p){
    for(jj in ii:p){ # symmetric hessian matrix
      hess.val <- -crossprod(d2mean[, ii, jj], SigX)
      hess.val <- hess.val + crossprod(SigMu[,ii], Sig2X[, jj])
      hess.val <- hess.val - crossprod(SigMu[,jj], Sig2X[, ii])
      hess.val <- hess.val + crossprod(dmean[, ii], Toep$Solve(dmean[, jj]))
      hess.val <- hess.val + crossprod(Sig2X[, jj], Toep$Aolve(Sig2X[, ii]))
      hess.val <- hess.val - Toep$TraceProd(d2acf[, ii, jj]) / 2
      hess.val <- hess.val + Toep$TraceDeriv(dacf[, jj], dacf[, ii]) / 2
      hess[ii, jj] <- hess.val
    }
  }
  hess <- hess + Sigd2X / 2
  hess
}
