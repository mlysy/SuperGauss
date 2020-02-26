#' Constructor and methods for `NormalToeplitz` objects.
#'
#' @details It is assumed that the autocorrelation of the [`Toeplitz`] object defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values typically contain `NaN`s).
#' @include Toeplitz.R
#' @export
.NormalToeplitz <- setRefClass("NormalToeplitz",
                               fields = list(cpp_ptr = "externalptr",
                                             size = "numeric",
                                             npara = "numeric"))
.NormalToeplitz$lock("cpp_ptr") # locked fields
.NormalToeplitz$lock("size")
.NormalToeplitz$lock("npara")
# internal constructor
.NormalToeplitz$methods(initialize = function(n, p) {
  cpp_ptr <<- .NormalToeplitz_constructor(n, p)
  size <<- n
  npara <<- p
})

# exported constructor
#' @rdname NormalToeplitz-class
#' @param n Size of the Toeplitz matrix.
#' @param p Number of unknown parameters.
#' @return A `NormalToeplitz` object.
#' @export
NormalToeplitz <- function(n, p) {
  Nt <- .NormalToeplitz$new(n, p)
  Nt
}


#--- custom methods ------------------------------------------------------------

# Log-density function
.NormalToeplitz$methods(logdens = function(z, acf) {
  if(length(z) != size) {
    stop("z has wrong length.")
  }
  if(length(acf) != size) {
    stop("acf has wrong length.")
  }
  .NormalToeplitz_logdens(cpp_ptr, z, acf)
})

# Gradient of log-density with respect to parameters
.NormalToeplitz$methods(grad = function(z, dz, acf, dacf) {
  if(length(z) != size) {
    stop("z has wrong length.")
  }
  if(!all(dim(dz) == c(size, npara))) {
    stop("dz has wrong dimensions.")
  }
  if(length(acf) != size) {
    stop("acf has wrong length.")
  }
  if(!all(dim(dacf) == c(size, npara))) {
    stop("dacf has wrong dimensions.")
  }
  .NormalToeplitz_grad(cpp_ptr, z, dz, acf, dacf)
})

# Hessian of log-density with respect to parameters
.NormalToeplitz$methods(hess = function(z, dz, d2z, acf, dacf, d2acf) {

  if(length(z) != size) {
    stop("z has wrong length.")
  }
  if(!all(dim(dz) == c(size, npara))) {
    stop("dz has wrong dimensions.")
  }
  if(!all(dim(d2z) == c(size, npara, npara))) {
    stop("d2z has wrong dimensions.")
  }
  if(length(acf) != size) {
    stop("acf has wrong length.")
  }
  if(!all(dim(dacf) == c(size, npara))) {
    stop("dacf has wrong dimensions.")
  }
  if(!all(dim(d2acf) == c(size, npara, npara))) {
    stop("d2acf has wrong dimensions.")
  }
  .NormalToeplitz_hess(cpp_ptr, z, dz, d2z, acf, dacf, d2acf)
})


# Full gradient of log-density
.NormalToeplitz$methods(grad_full = function(z, acf) {

  if(length(z) != size) {
    stop("z has wrong length.")
  }
  if(length(acf) != size) {
    stop("acf has wrong length.")
  }

  .NormalToeplitz_grad_full(cpp_ptr, z, acf)
})

# now make methods display arguments
.DollarNames.NormalToeplitz <- function(x, pattern)
  grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
