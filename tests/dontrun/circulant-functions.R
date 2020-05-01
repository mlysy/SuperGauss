#--- Gaussian likelihoods with Circulant variance matrices ---------------------

require(R6)
require(fftw)
# columnwise fft and inverse fft for real inputs
fft <- function(x) {
  eV <- is.vector(x)
  y <- apply(as.matrix(x), 2, fftw::FFT)
  if(eV) y <- drop(y)
  y
}
ifft <- function(x) {
  eV <- is.vector(x)
  y <- apply(as.matrix(x), 2, fftw::FFT, inverse = TRUE)
  y <- Re(y)/nrow(y)
  if(eV) y <- drop(y)
  y
}

t2acf <- function(N, tacf) {
  n <- length(tacf)
  if(n != floor(N/2) + 1) stop("tacf has wrong length.")
  acf <- rep(NA, N)
  acf[1:n] <- tacf
  eN <- (2*n) == (N+2)
  id <- n - eN + (2:n-1)
  acf[n - eN + (2:n-1)] <- tacf[n:2]
  acf
}

Circulant <- R6Class(
  classname = "Circulant",
  private = list(
    N_ = NA,
    n_ = NA,
    acf_ = NULL,
    psd_ = NULL,
    ldet_ = NULL,
    #' @description Precomputations for power product.
    prod_setup = function() {
      if(!self$has_acf()) stop("acf has not been set.")
      private$psd_ <- Re(fft(private$acf_))
    },
    #' @description FFT multiplication.
    #'
    #' Calculates `z = y %*% x`, where `yf = fft(y)`.
    #' Return value `z` has same dimensions as `x`.
    fft_mult = function(x, yf) {
      vx <- is.vector(x)
      if(vx) x <- as.matrix(x)
      if(!is.matrix(x) || (nrow(x) != private$N_)) {
        stop("x has incorrect dimensions.")
      }
      ## xf <- apply(x, 2, fftw::FFT)
      xf <- fft(x)
      zf <- yf * xf
      ## z <- Re(apply(zf, 2, fftw::FFT, inverse = TRUE))/private$N_
      z <- ifft(zf)
      if(vx) z <- drop(z)
      z
    }
  ),
  public = list(
    #' @description Constructor
    #'
    #' @param N Size of Circulant matrix.
    #' @param tacf Truncated autocorrelation vector of length `n = floor(N/2) + 1`.
    #'
    #' @return A `Circulant` object.
    initialize = function(N, tacf) {
      private$N_ <- N
      private$n_ <- floor(N/2) + 1
      if(!missing(tacf)) self$set_acf(tacf)
    },
    #' @description Set the autocorrelation of the Circulant matrix.
    #'
    #' @param tacf Truncated autocorrelation vector of length `n = floor(N/2) + 1`.
    set_acf = function(tacf) {
      if(length(tacf) != private$n_) {
        stop("tacf has wrong length.")
      }
      private$acf_[1:private$n_] <- tacf
      if(private$N_ > 1) {
        eN <- (2*private$n_) == (private$N_+2)
        id <- private$n_ - eN + (2:private$n_-1)
        private$acf_[id] <- tacf[private$n_:2]
      }
    },
    #' @description Get the autocorrelation of the Circulant matrix.
    #'
    #' @return The autocorrelation vector of length `N`.
    get_acf = function() {
      private$acf_
    },
    #' @description Check if the autocorrelation has been set.
    #'
    has_acf = function() {
      !is.null(private$acf_)
    },
    #' @description Get the truncated autocorrelation of the Circulant matrix.
    #'
    #' @return The truncated autocorrelation vector of length `n = floor(N/2) + 1`.
    get_tacf = function() {
      private$acf_[1:private$n_]
    },
    #' @description Power product.
    #'
    #' Calculates `Ct^pow %*% x`.
    #'
    #' @param x Vector or matrix.  If missing defaults to `diag(N)`.
    #' @param pow Power parameter (scalar).
    pow_prod = function(x, pow = 1) {
      if(is.null(private$psd_)) private$prod_setup()
      if(missing(x)) x <- diag(private$N_)
      yf <- private$psd_^pow
      private$fft_mult(x, yf)
    },
    #' @description Solve circulant system.
    #'
    #' Calculates `Ct^{-1} x`.
    #' @param x Vector or matrix.  If missing defaults to `diag(N)`.
    solve = function(x) {
      if(is.null(private$psd_)) private$prod_setup()
      if(missing(x)) x <- diag(private$N_)
      private$fft_mult(x, yf = 1/private$psd_)
    },
    #' @description Product with circulant matrix.
    #'
    #' Calculates `Ct %*% x`.
    #' @param x Vector or matrix.  If missing defaults to `diag(N)`.
    prod = function(x) {
      if(is.null(private$psd_)) private$prod_setup()
      if(missing(x)) x <- diag(private$N_)
      private$fft_mult(x, yf = private$psd_)
    },
    #' @description Log-determinant.
    #'
    #' Calculates `log|Ct|`.
    #'
    log_det = function() {
      if(is.null(private$ldet_)) {
        if(is.null(private$psd_)) private$prod_setup()
        private$ldet_ <- sum(log(private$psd_))
      }
      private$ldet_
    }
  )
)


N <- 6
n <- floor(N/2) + 1
Ct <- Circulant$new(N = N, tacf = 1:n)
Ct$get_acf()

#--- power products with circulants --------------------------------------------


N <- 6
n <- floor(N/2) + 1
Ct <- Circulant$new(N = N, tacf = n:1)
Ct$get_acf()

N <- sample(5:20,1)
p <- sample(1:5, 1)
n <- floor(N/2) + 1
tacf <- rnorm(n)
Ct <- Circulant$new(N = N, tacf = tacf)
x <- matrix(rnorm(N*p), N, p)
y <- toeplitz(Ct$get_acf()) %*% x
y2 <- Ct$pow_prod(x)
range(y-y2)
y <- solve(toeplitz(Ct$get_acf()), x)
y2 <- Ct$pow_prod(x, pow = -1)
range(y - y2)

#--- circulant convolution -----------------------------------------------------

require(mvtnorm)
require(numDeriv)

N <- 5
x <- rnorm(N)

(1:N-1 + 0) %% N + 1

#' Circulant matrix with first row `x`.
circulant <- function(x) {
  N <- length(x)
  t(sapply(1:N-1, function(ii) x[(1:N-1 - ii) %% N + 1]))
}

#' Shift matrix.
shift <- function(N, i) {
  circulant(as.numeric(1:N == i))
}

#' Calculate the vector `(x' S_0 y, x' S_1 y, ..., x' S_N y)`, where `S_i` is the `i`th circulant shift matrix.
shift_mult <- function(x, y) {
  rev(ifft(fft(x) * fft(rev(y))))
  ## N <- length(x)
  ## rev(Re(FFT(FFT(x) * FFT(rev(y)), inverse = TRUE)/N))
}

shift(5, 2)

N <- 10
x <- rnorm(N)
y <- rnorm(N)

c1 <- sapply(1:N, function(i) crossprod(x, shift(N, i) %*% y))
c2 <- rev(convolve(x, rev(y), conj = FALSE))
c3 <- shift_mult(x, y)
range(c1 - c2)
range(c1 - c3)

# ok let's try the gradient

circ_ldens <- function(z, nu) {
  Ct <- Circulant$new(N = length(z), tacf = nu)
  ld <- sum(z*Ct$solve(z)) + Ct$log_det()
  -.5 * ld
}

N <- 15
n <- floor(N/2)+1
psd <- rexp(n)
gam <- Re(FFT(t2acf(N, psd)))
nu <- gam[1:n]
z <- rnorm(N)
circ_ldens(z, nu) - dmvnorm(z, sigma = toeplitz(t2acf(N, nu)), log = TRUE)

circ_ip <- function(z, nu) {
  Ct <- Circulant$new(N = length(z), tacf = nu)
  sum(z*Ct$solve(z))
}

circ_tr <- function(N, nu) {
  Ct <- Circulant$new(N = N, tacf = nu)
  Ct$log_det()
}

N <- sample(2:20,1)
n <- floor(N/2)+1
## nu <- Re(FFT(t2acf(N, rexp(n)*N), inverse = TRUE)[1:n]/N)
psd <- rexp(n)*N
nu <- ifft(t2acf(N, psd))[1:n]
z <- rnorm(N)
g1 <- grad(circ_ip, x = nu, z = z)
y <- solve(toeplitz(t2acf(N, nu)), z)
g2 <- -shift_mult(y, y)
if(N > 2) {
  eN <- (2*n) == (N+2)
  id <- 2:(n-eN)
  g2[id] <- 2*g2[id]
}
range(g1-g2[1:n])


N <- sample(5:15, 1)
n <- floor(N/2)+1
psd <- rexp(n)*N
## gam <- Re(FFT(t2acf(N, psd), inverse = TRUE))/N
gam <- ifft(t2acf(N, psd))
nu <- gam[1:n]
V1 <- solve(toeplitz(gam))
V2 <- toeplitz(ifft(t2acf(N, 1/psd)))
range(V1 - V2)

N <- sample(5:20,1)
n <- floor(N/2)+1
psd <- rexp(n)*N
gam <- ifft(t2acf(N, psd))
nu <- gam[1:n]
g1 <- grad(circ_tr, x = nu, N = N)
## g2 <- Re(FFT(1/Re(FFT(t2acf(N, nu))), inverse = TRUE))
g2 <- ifft(1/fft(t2acf(N, nu))) * N
if(N > 2) {
  eN <- (2*n) == (N+2)
  id <- 2:(n-eN)
  g2[id] <- 2*g2[id]
}
range(g1-g2[1:n])
