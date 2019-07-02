#' @title Constructor and methods for PCG objects.
#'
#' @description The \code{PCG} class contains efficient methods for linear algebra with symmetric positive definite (i.e., variance) Toeplitz matrices.
#'
#' @section Methods:
#' If \code{P1} is a \code{PCG} object, then:
#' \describe{
#' \item{\code{solve(Toep, X)}, \code{solve(Toep)}}{Solves Toeplitz systems of equations.  When second argument is missing, returns the inverse of the Toeplitz matrix.}
#' }
#' @details It is assumed that the autocorrelation of the \code{Toeplitz} object defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values typically contain \code{NaN}s).
#' @examples
#' # construction
#' N <- 5
#' acf <- exp(-(1:N))
#' P1 <- PCG(n = N)
#' 
#' # linear algebra
#' X <- rnorm(N)
#' P1$solve(acf = acf, y = X)
#' @export
.PCG <- setRefClass("PCG",
                         fields = list(cpp_ptr = "externalptr",
                                       size = "numeric"))
.PCG$lock("cpp_ptr") # locked fields
.PCG$lock("size")
# internal constructor
.PCG$methods(initialize = function(n) {
  cpp_ptr <<- .PCG_constructor(n)
  size <<- n
})

# exported constructor
#' @rdname PCG-class
#' @param n Size of the Toeplitz matrix.
#' @return A \code{PCG} object.
#' @export
PCG <- function(n) {
  P1 <- .PCG$new(n)
  P1
}

#--- custom methods ------------------------------------------------------------

# Solve
.PCG$methods(solve = function(acf, y, tol) {
  .PCG_solve(cpp_ptr, acf, y, tol)
})
