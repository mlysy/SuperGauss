check_acf <- function(acf, N) {
  if(!is.numeric(acf) || !is.vector(acf)) {
    stop("acf must be a numeric vector.")
  }
  if(length(acf) != N) {
    stop("acf and X have incompatible dimensions.")
  }
}

# extract number of parameters from these values
# returns an error if this is not possible or
# there are conflicting number of parameters
.get.p <- function(dmu, dacf) {
  p <- rep(NA, 2)
  if(!missing(dmu)) {
    p[1] <- if(is.vector(dmu)) length(dmu) else ncol(dmu)
  }
  if(!missing(dacf)) {
    if(!is.matrix(dacf)) {
      stop("dacf must be a matrix.")
    } else {
      p[2] <- ncol(dacf)
    }
  }
  if(all(is.na(p))) {
    stop("Must provide either dmu or dacf.")
  }
  # check this elsewhere
  ## if(!anyNA(p)) {
  ##   if(p[1] != p[2]) {
  ##     stop("dmu and dacf have incompatible dimensions.")
  ##   }
  ## }
  p[!is.na(p)][1]
}

# mu can be missing, a scalar, or a length-N vector
# dmu can be missing, a length-p vector, or an N x p matrix
# d2mu can be missing, a p x p matrix, or an N x p x p array
# returns an error if mu is missing but dmu or d2mu are not
.format.mu <- function(mu, dmu, d2mu, N, p, grad.only = TRUE) {
  if(!missing(dmu) && missing(mu)) {
    stop("Cannot provide dmu without mu.")
  }
  if(!missing(d2mu) && (missing(mu) || missing(dmu))) {
    stop("Cannot provide d2mu without both mu and dmu.")
  }
  if(missing(mu)) {
    mu <- 0
  }
  # format mu
  if(length(mu) == 1) mu <- rep(mu, N)
  # format dmu
  if(missing(dmu)) dmu <- matrix(0, N, p)
  if(is.vector(dmu)) dmu <- matrix(dmu, N, length(dmu), byrow = TRUE)
  if(nrow(dmu) != N) {
    stop("dmu and X have incompatible dimensions.")
  }
  if(ncol(dmu) != p) {
    stop("dmu and dacf have incompatible dimensions.")
  }
  out <- list(mu = mu, dmu = dmu)
  # format d2mu
  if(!grad.only) {
    if(missing(d2mu)) d2mu <- array(0, dim = c(N, p, p))
    if(is.matrix(d2mu)) {
      d2mu <- array(c(d2mu), dim = c(nrow(d2mu), ncol(d2mu), N))
      d2mu <- aperm(d2mu, perm = c(3, 1, 2))
    }
    if(dim(d2mu)[1] != N) {
      stop("d2mu and X have incompatible dimensions.")
    }
    if(!identical(dim(d2mu)[2:3], c(p,p))) {
      stop("d2mu and dmu/dacf have incompatible dimensions.")
    }
    out <- c(out, list(d2mu = d2mu))
  }
  out
}

.format.dacf <- function(dacf, d2acf, N, p, grad.only = TRUE) {
  if(!missing(d2acf) && missing(dacf)) {
    stop("Cannot provide d2acf without dacf.")
  }
  if(missing(dacf)) dacf <- matrix(0, N, p)
  if(is.vector(dacf)) dacf <- matrix(dacf, N, length(dacf), byrow = TRUE)
  if(nrow(dacf) != N) stop("dacf and X have incompatible dimensions.")
  if(ncol(dacf) != p) stop("dmu and dacf have incompatible dimensions.")
  out <- list(dacf = dacf)
  if(!grad.only) {
    if(missing(d2acf)) d2acf <- array(0, dim = c(N, p, p))
    if(is.matrix(d2acf)) {
      d2acf <- array(c(d2acf), dim = c(nrow(d2acf), ncol(d2acf), N))
      d2acf <- aperm(d2acf, perm = c(3, 1, 2))
    }
    if(dim(d2acf)[1] != N) stop("d2acf and X have incompatible dimensions.")
    if(!identical(dim(d2acf)[2:3], c(p,p))) {
      stop("d2acf and dmu/dacf have incompatible dimensions.")
    }
    out <- c(out, list(d2acf = d2acf))
  }
  out
}

# acf should be Toeplitz object of size N
.format.acf <- function(acf, N) {
  if(is.Toeplitz(acf)) {
    # is Toeplitz
    if(acf$size() != N) {
      stop("acf and X have incompatible dimensions.")
    }
  } else if(is.vector(acf)) {
    if(length(acf) != N) {
      stop("acf and X have incompatible dimensions.")
    } else {
      acf2 <- acf
      acf <- Toeplitz$new(N)
      acf$set_acf(acf2)
    }
  } else {
    stop("acf should be either a vector or a Toeplitz object.")
  }
  acf
}
