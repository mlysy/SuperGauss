#' Constructor and methods for `NormalToeplitz` objects.
#'
#' @details It is assumed that the autocorrelation of the [`Toeplitz`] object defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values typically contain `NaN`s).
#' @include Toeplitz.R
#' @export
.NormalToeplitz <- setRefClass("NormalToeplitz",
                               fields = list(cpp_ptr = "externalptr",
                                             size = "numeric"))
.NormalToeplitz$lock("cpp_ptr") # locked fields
.NormalToeplitz$lock("size")
## .NormalToeplitz$lock("npara")
# internal constructor
.NormalToeplitz$methods(initialize = function(n) {
  cpp_ptr <<- .NormalToeplitz_constructor(n)
  size <<- n
  ## npara <<- p
})

# exported constructor
#' @rdname NormalToeplitz-class
#' @param n Size of the Toeplitz matrix.
#' @return A `NormalToeplitz` object.
#' @export
NormalToeplitz <- function(n) {
  Nt <- .NormalToeplitz$new(n)
  Nt
}


#--- custom methods ------------------------------------------------------------

# Log-density function.
.NormalToeplitz$methods(logdens = function(z, acf) {
  check_var(z, N = size, varname = "z")
  check_var(acf, N = size, varname = "acf")
  .NormalToeplitz_logdens(cpp_ptr, z, acf)
})

# Gradient of log-density with respect to parameters
.NormalToeplitz$methods(grad = function(z, dz, acf, dacf) {
  check_var(z, dx = dz, N = size, varname = "z")
  check_var(acf, dx = dacf, N = size, varname = "acf")
  if(ncol(dz) != ncol(dacf)) {
    stop("dz and dacf must have the same number of columns.")
  }
  ntheta <- ncol(dz)
  .NormalToeplitz_grad(cpp_ptr, z, dz, acf, dacf, ntheta)
})

# Hessian of log-density with respect to parameters
.NormalToeplitz$methods(hess = function(z, dz, d2z, acf, dacf, d2acf) {
  check_var(z, dx = dz, d2x = d2z, N = size, varname = "z")
  check_var(acf, dx = dacf, d2x = d2acf, N = size, varname = "acf")
  if(ncol(dz) != ncol(dacf)) {
    stop("dz and dacf must have the same number of columns.")
  }
  if(!all(dim(d2z) == dim(d2acf))) {
    stop("d2z and d2acf must have the same dimensions.")
  }
  ntheta <- ncol(dz)
  .NormalToeplitz_hess(cpp_ptr, z, dz,
                       matrix(d2z, size, ntheta*ntheta),
                       acf, dacf,
                       matrix(d2acf, size, ntheta*ntheta), ntheta)
})


# Full gradient of log-density
.NormalToeplitz$methods(grad_full = function(z, acf) {
  check_var(z, N = size, varname = "z")
  check_var(acf, N = size, varname = "acf")
  .NormalToeplitz_grad_full(cpp_ptr, z, acf)
})

# now make methods display arguments
.DollarNames.NormalToeplitz <- function(x, pattern) {
  grep(pattern, getRefClass(class(x))$methods(), value = TRUE)
}

#--- formatting functions ------------------------------------------------------

# Check inputs to NormalToeplitz
#
# @param x Main input variable.
# @param dx Optional gradient variable.
# @param d2x Optional Hessian variable.
# @param varname Name of variable.
#
# @details If only `x` provided, checks length.  If `dx` provided, checks that its a matrix with correct dimensions.  If `d2x` is provided, check that its' an array with dimensions compatible with `x` and `dx`.
check_var <- function(x, dx, d2x, N, varname) {
  if(!(is.vector(x) && is.numeric(x))) {
    stop(paste0(varname, " must be a numeric vector."))
  }
  if(length(x) != N) {
    stop(paste0(varname, "has wrong length."))
  }
  if(!missing(dx)) {
    if(!(is.matrix(dx) && is.numeric(dx))) {
      stop(paste0("d", varname, " must be a numeric matrix."))
    }
    if(nrow(dx) != N) {
      stop(paste0("d", varname, " has wrong number of rows."))
    }
    p <- ncol(dx)
    if(!missing(d2x)) {
      if(!(is.array(d2x) && (length(dim(d2x)) == 3) && is.numeric(d2x))) {
        stop(paste0("d2", varname, " must be a numeric 3D array."))
      }
      if(!all(dim(d2x) == c(N, p, p))) {
        stop(paste0("d2", varname, " has incompatible dimensions with d",
                    varname, " and ", varname, "."))
      }
    }
    if(!missing(d2x) && missing(dx)) {
      stop(paste0("Cannot provide d2", varname, " without d", varname, "."))
    }
  }
}
