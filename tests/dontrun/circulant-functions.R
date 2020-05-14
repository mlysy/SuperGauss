#--- Gaussian likelihoods with Circulant variance matrices ---------------------

require(R6)
require(fftw)

# columnwise fft and inverse fft for real inputs
fft <- function(x) {
  eV <- is.vector(x)
  x <- as.matrix(x)
  if(nrow(x) == 1) {
    y <- x + 0i
  }  else {
    y <- apply(as.matrix(x), 2, fftw::FFT)
  }
  if(eV) y <- drop(y)
  y
}
ifft <- function(x) {
  eV <- is.vector(x)
  x <- as.matrix(x)
  if(nrow(x) == 1) {
    y <- x + 0i
  }  else {
    y <- apply(as.matrix(x), 2, fftw::FFT, inverse = TRUE)
  }
  y <- Re(y)/nrow(y)
  if(eV) y <- drop(y)
  y
}

# convert a truncated acf to full acf.
t2acf <- function(N, tacf) {
  n <- length(tacf)
  if(n != floor(N/2) + 1) stop("tacf has wrong length.")
  acf <- rep(NA, N)
  acf[1:n] <- tacf
  if(N > 1) {
    eN <- (2*n) == (N+2)
    id <- n - eN + (2:n-1)
    acf[n - eN + (2:n-1)] <- tacf[n:2]
  }
  acf
}

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

# circulant matrix class.
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

# Normal with circulant variance matrix.
NormalCirculant <- R6Class(
  classname = "NormalCirculant",
  private = list(
    Ct_ = NULL,
    N_ = NA,
    n_ = NA
  ),
  public = list(
    initialize = function(N) {
      private$N_ <- N
      private$n_ <- (N %/% 2) + 1
      private$Ct_ <- Circulant$new(N = N)
    },
    logdens = function(z, tacf) {
      private$Ct_$set_acf(tacf)
      ld <- sum(z*private$Ct_$solve(z)) + private$Ct_$log_det()
      -.5 * (ld + private$N_ * log(2*pi))
    },
    grad_full = function(z, tacf) {
      private$Ct_$set_acf(tacf)
      # (negative) gradient wrt z
      dz <- private$Ct_$solve(z)
      # gradient wrt tacf
      dip <- shift_mult(dz, dz) # gradient of inner product
      dtr <- ifft(1/fft(t2acf(private$N_, tacf))) * private$N_ # gradient of trace
      if(private$N_ > 2) {
        eN <- (2*private$n_) == (private$N_+2)
        id <- 2:(private$n_-eN)
        dip[id] <- 2*dip[id]
        dtr[id] <- 2*dtr[id]
      }
      dtacf <- .5 * (dip[1:private$n_] - dtr[1:private$n_])
      list(dz = -dz, dtacf = dtacf)
    }
  )
)

