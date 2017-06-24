#' Constructor and methods for Toeplitz matrix objects.
#'
#' The \code{Toeplitz_Matrix} class contains efficient methods for linear algebra with symmetric positive definite (i.e., variance) Toeplitz matrices.
#'
#' @aliases setAcf getAcf traceT2 traceT4 show.Toeplitz_Matrix %*% determinant solve
#' @section Methods:
#' If \code{Toep} is a \code{Toeplitz_Matrix} object with first row/column given by \code{acf}, then:
#' \describe{
#' \item{\code{Toep$setAcf(acf)}}{Sets the autocorrelation of the matrix.}
#' \item{\code{Toep$getAcf()}}{Gets the autocorrelation of the matrix.}
#' \item{\code{nrow(Toep)}, \code{ncol(Toep)}, \code{dim(Toep)}}{Selected dimension(s) of the matrix.}
#' \item{\code{Toep \%*\% X}, \code{X \%*\% Toep}}{Toeplitz-Matrix and Matrix-Toeplitz multiplication.  Also works if \code{X} is a vector.}
#' \item{\code{solve(Toep, X)}, \code{solve(Toep)}}{Solves Toeplitz systems of equations.  When second argument is missing, returns the inverse of the Toeplitz matrix.}
#' \item{\code{determinant(Toep)}}{Log-determinant of the Toeplitz matrix, i.e., same thing as \code{determinant(toeplitz(acf))$modulus}.}
#' \item{\code{Toep$traceT2(acf2)}}{Computes the trace of \code{solve(toeplitz(acf), toeplitz(acf2))}.  This is used in the computation of the gradient of Gaussian likelihoods with Toeplitz variance matrix.}
#' \item{\code{Toep$traceT4(acf2, acf3)}}{Computes the trace of \code{solve(toeplitz(acf), toeplitz(acf2)) \%*\% solve(toeplitz(acf), toeplitz(acf3))}.  This is used in the computation of the Hessian of Gaussian likelihoods with Toeplitz variance matrix.}
#' }
#' @details It is assumed that the autocorrelation of the \code{Toeplitz_Matrix} defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values are typically \code{NaN}s).
#' @examples
#' # construction
#' acf <- exp(-(1:5))
#' Toep <- Toeplitz(acf = acf)
#' \dontrun{
#' # alternatively, can allocate space first
#' Toep <- Toeplitz(n = length(acf))
#' Toep$setAcf(acf = acf)
#' }
#' dim(Toep)
#' Toep # show method
#' Toep$getAcf() # extract the acf
#'
#' # linear algebra
#' X <- matrix(rnorm(10), 5, 2)
#' Toep %*% X
#' t(X) %*% Toep
#' solve(Toep, X)
#' determinant(Toep) # log-determinant
#' @export
.Toeplitz <- setRefClass("Toeplitz_Matrix",
                         fields = list(cpp_ptr = "externalptr",
                                       size = "numeric"))
.Toeplitz$lock("cpp_ptr") # locked fields
.Toeplitz$lock("size")
# internal constructor
.Toeplitz$methods(initialize = function(n) {
  cpp_ptr <<- .Toeplitz_constructor(n)
  size <<- n
})

# exported constructor
#' @rdname Toeplitz_Matrix-class
#' @param n size of the Toeplitz matrix.
#' @param acf autocorrelation vector of Toeplitz matrix.
#' @return A \code{Toeplitz_Matrix} object.
#' @export
Toeplitz <- function(n, acf) {
  if(missing(n)){
    n <- length(acf)
  }
  Tz <- .Toeplitz$new(n)
  if(!missing(acf)) {
    Tz$setAcf(acf)
  }
  Tz
}


#--- custom methods ------------------------------------------------------------

# setter
.Toeplitz$methods(setAcf = function(acf) {
  if(length(acf) != size) {
    stop("acf has wrong length.")
  }
  .Toeplitz_setAcf(cpp_ptr, acf)
})

# getter
.Toeplitz$methods(getAcf = function() {
  .Toeplitz_getAcf(cpp_ptr)
})

# traceT2
.Toeplitz$methods(traceT2 = function(acf2) {
  if(!.Toeplitz_hasAcf(cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(length(acf2) != size) {
    stop("acf2 has wrong length.")
  }
  .Toeplitz_traceT2(cpp_ptr, acf2)
})

# traceT4
.Toeplitz$methods(traceT4 = function(acf2, acf3) {
  if(!.Toeplitz_hasAcf(cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(length(acf2) != size) {
    stop("acf2 has wrong length.")
  }
  if(length(acf3) != size) {
    stop("acf3 has wrong length.")
  }
  .Toeplitz_traceT4(cpp_ptr, acf2, acf3)
})

#--- generic methods -----------------------------------------------------------

# show
#' @export
setMethod("show", "Toeplitz_Matrix", function(object) {
  if(.Toeplitz_hasAcf(object$cpp_ptr)) {
    obj.acf <- object$getAcf()[1:min(6, object$size)]
    obj.acf <- signif(obj.acf, digits = 3)
    if(object$size > 6) obj.acf <- c(obj.acf, "...")
  } else {
    obj.acf <- "NULL"
  }
  cat("Toeplitz_Matrix of size", object$size, "\n",
      "acf: ", obj.acf, "\n")
})

# ncol
#' @export
setMethod("ncol", "Toeplitz_Matrix", function(x){
  x$size
})

# nrow
#' @export
setMethod("nrow", "Toeplitz_Matrix", function(x){
  x$size
})

# dim
#' @export
setMethod("dim", "Toeplitz_Matrix", function(x){
  rep(x$size, 2)
})

# Toeplitz-Matrix multiplication
#' @export
setMethod("%*%", signature(x = "Toeplitz_Matrix", y = "ANY"), function(x, y) {
  if(!.Toeplitz_hasAcf(x$cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(is.vector(y)) y <- as.matrix(y)
  if(!is.matrix(y)) {
    stop("Second argument should be a matrix or vector.")
  }
  if(nrow(y) != x$size) {
    stop("Toeplitz_Matrix and second argument are non-conformable.")
  }
  .Toeplitz_Multiply(x$cpp_ptr, y)
})

# Matrix-Toeplitz multiplication
setMethod("%*%", signature(x = "ANY", y = "Toeplitz_Matrix"), function(x, y) {
  if(!.Toeplitz_hasAcf(y$cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(is.vector(x)) x <- as.matrix(x)
  if(!is.matrix(x)) {
    stop("First argument should be a matrix or vector.")
  }
  if(ncol(x) != y$size) {
    stop("First argument and Toeplitz_Matrix are non-conformable.")
  }
  t(.Toeplitz_Multiply(y$cpp_ptr, t(x)))
})

# determinant
#' @export
setMethod("determinant", "Toeplitz_Matrix",
          function(x, logarithm = TRUE, ...) {
  if(!.Toeplitz_hasAcf(x$cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  ldT <- .Toeplitz_Determinant(x$cpp_ptr)
  if(!logarithm) {
    ldT <- exp(ldT)
  }
  ldT
})

# solve
#' @export
setMethod("solve", "Toeplitz_Matrix", function(a, b, ...) {
  if(!.Toeplitz_hasAcf(a$cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(missing(b)) b <- diag(a$size)
  if(is.vector(b)) b <- as.matrix(b)
  if(!is.matrix(b)) {
    stop("b must be a matrix or vector.")
  }
  if(nrow(b) != a$size) {
    stop("a and b are non-conformable.")
  }
  .Toeplitz_Solve(a$cpp_ptr, b)
})


# new Toeplitz class.
# write doc for this later, but basic idea is that class is implemented with XPtr instead of Modules.  The reason for this is that ultimately there is overhead (since Toeplitz class in C++ needs to be wrapped manually anyways), and also the exposed C++ members are highly vulnerable to incorrect inputs.

# so here are the list of mewthods in the R Toeplitz object:
# * constructor.  C++ destructor called automatically when R object is deleted
# * setAcf/getAcf.  Getter and setter.  The former can now return an error if it has the wrong size.
# * nrow/ncol/dim.  These are generic methods.
# * solve.  generic method.
# * determinant/det.  generic method.
# * mult.  generic method.  only works one way though...
# * traceT2/traceT4.  these are both non-generics.

# ok let's define some of these.
