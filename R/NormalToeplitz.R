#' @title Constructor and methods for NormalToeplitz objects.
#'
#' @description 
#'
#' @details It is assumed that the autocorrelation of the \code{Toeplitz} object defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values typically contain \code{NaN}s).
#' @examples
#' @export
.NormalToeplitz <- setRefClass("NormalToeplitz",
                               fields = list(cpp_ptr = "externalptr",
                                             size = "numeric", 
                                             npara = "numeric"))
.NormalToeplitz$lock("cpp_ptr") # locked fields
.NormalToeplitz$lock("size")
.NormalToeplitz$lock("npara")
# internal constructor
.Toeplitz$methods(initialize = function(n, p) {
  cpp_ptr <<- .Toeplitz_constructor(n, p)
  size <<- n
  npara <<- p
})

# exported constructor
#' @rdname NormalToeplitz-class
#' @param n Size of the Toeplitz matrix.
#' @param p Number of unknown parameters.
#' @param hasToep Flag indicating whether Toeplitz member is pre-generated.
#' @return A \code{NormalToeplitz} object.
#' @export
NormalToeplitz <- function(n, p) {
  Nt <- .NormalToeplitz$new(n, p)
  Nt
}


#--- generic methods -----------------------------------------------------------

# Log-density function
#' @export
setMethod("logdens", signature(x = "NormalToeplitz", y = "ANY", z = "ANY"), 
          function(x, y, z) {
  .NormalToeplitz_logdens(x$cpp_ptr, y, z)
})

# Gradient of log-density with respect to parameters
#' @export
setMethod("grad", signature(x = "NormalToeplitz", z = "ANY", dzdt = "ANY", 
                            y = "ANY", dacfdt = "ANY"), 
          function(x, z, dzdt, y, dacfdt) {
  .NormalToeplitz_grad(x$cpp_ptr, z, dzdt, y, dacfdt)
})

# Hessian of log-density with respect to parameters
#' @export
setMethod("hess", signature(x = "NormalToeplitz", z = "ANY", dzdt = "ANY", d2zdt = "ANY",
                            y = "ANY", dacfdt = "ANY", d2acfdt = "ANY"), 
          function(x, z, dzdt, d2zdt, y, dacfdt, d2acfdt) {
  .NormalToeplitz_hess(x$cpp_ptr, z, dzdt, d2zdt, y, dacfdt, d2acfdt)
})

# Full gradient of log-density
#' @export
setMethod("grad_full", signature(x = "NormalToeplitz", z = "ANY", y = "ANY"), 
          function(x, z, y) {
  .NormalToeplitz_grad_full(x$cpp_ptr, z, y)
})

#' # Overloaded Log-density function
#' #' @export
#' setMethod("logdens", signature(x = "NormalToeplitz", y = "ANY", z = "Toeplitz"), 
#'           function(x, y, z) {
#'             .NormalToeplitz_logdens_(x$cpp_ptr, y, z$cpp_ptr)
#'           })
#' 
#' # Overloaded gradient of log-density with respect to parameters
#' #' @export
#' setMethod("grad", signature(x = "NormalToeplitz", z = "ANY", dzdt = "ANY", 
#'                             y = "Toeplitz", dacfdt = "ANY"), 
#'           function(x, z, dzdt, y, dacfdt) {
#'             .NormalToeplitz_grad_(x$cpp_ptr, z, dzdt, y$cpp_ptr, dacfdt)
#'           })
#' # Overloaded Hessian of log-density with respect to parameters
#' #' @export
#' setMethod("hess", signature(x = "NormalToeplitz", z = "ANY", dzdt = "ANY", d2zdt = "ANY",
#'                             y = "Toeplitz", dacfdt = "ANY", d2acfdt = "ANY"), 
#'           function(x, z, dzdt, d2zdt, y, dacfdt, d2acfdt) {
#'             .NormalToeplitz_hess_(x$cpp_ptr, z, dzdt, d2zdt, y$cpp_ptr, dacfdt, d2acfdt)
#'           })
#' # Overloaded full gradient of log-density
#' #' @export
#' setMethod("grad_full", signature(x = "NormalToeplitz", z = "ANY", y = "Toeplitz"), 
#'           function(x, z, y) {
#'             .NormalToeplitz_grad_full_(x$cpp_ptr, z, y$cpp_ptr)
#'           })

# now make methods display arguments
.DollarNames.NormalToeplitz <- function(x, pattern)
  grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
