#' @title Constructor and methods for Toeplitz matrix objects.
#'
#' @description The \code{Toeplitz} class contains efficient methods for linear algebra with symmetric positive definite (i.e., variance) Toeplitz matrices.
#'
#' @aliases setAcf getAcf traceT2 traceT4 show.Toeplitz %*% determinant solve %*%,ANY,Toeplitz-method %*%,Toeplitz,ANY-method determinant,Toeplitz-method dim,Toeplitz-method ncol,Toeplitz-method nrow,Toeplitz-method show,Toeplitz-method solve,Toeplitz-method
#' @section Methods:
#' If \code{Toep} is a \code{Toeplitz} object with first row/column given by \code{acf}, then:
#' \describe{
#' \item{\code{Toep$setAcf(acf)}}{Sets the autocorrelation of the matrix.}
#' \item{\code{Toep$getAcf()}}{Gets the autocorrelation of the matrix.}
#' \item{\code{nrow(Toep)}, \code{ncol(Toep)}, \code{dim(Toep)}}{Selected dimension(s) of the matrix.}
#' \item{\code{Toep \%*\% X}, \code{X \%*\% Toep}}{Toeplitz-Matrix and Matrix-Toeplitz multiplication.  Also works if \code{X} is a vector.}
#' \item{\code{solve(Toep, X)}, \code{solve(Toep)}}{Solves Toeplitz systems of equations.  When second argument is missing, returns the inverse of the Toeplitz matrix.}
#' \item{\code{determinant(Toep)}}{Log-determinant of the Toeplitz matrix, i.e., same thing as \code{log(det(toeplitz(acf)))}.}
#' \item{\code{Toep$traceT2(acf2)}}{If \code{T1 == toeplitz(acf)} and \code{T2 == toeplitz(acf2)}, computes the trace of \code{solve(T1, T2)}.  This is used in the computation of the gradient of Gaussian likelihoods with Toeplitz variance matrix.}
#' \item{\code{Toep$traceT4(acf2, acf3)}}{If \code{T1 == toeplitz(acf)}, \code{T2 == toeplitz(acf2)}, and \code{T3 == toeplitz(acf3)}, computes the trace of \code{solve(T1, T2) \%*\% solve(T1, T3)}.  This is used in the computation of the Hessian of Gaussian likelihoods with Toeplitz variance matrix.}
#' }
#' @details It is assumed that the autocorrelation of the \code{Toeplitz} object defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values typically contain \code{NaN}s).
#' @examples
#' # construction
#' acf <- exp(-(1:5))
#' Toep <- Toeplitz(acf = acf)
#' # alternatively, can allocate space first
#' Toep <- Toeplitz(n = length(acf))
#' Toep$setAcf(acf = acf)
#'
#' dim(Toep) # == c(nrow(Toep), ncol(Toep))
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
.Toeplitz <- setRefClass("Toeplitz",
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
#' @rdname Toeplitz-class
#' @param n Size of the Toeplitz matrix.
#' @param acf Autocorrelation vector of Toeplitz matrix.
#' @return A \code{Toeplitz} object.
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

#--- generic methods -----------------------------------------------------------

# show
#' @export
setMethod("show", "Toeplitz", function(object) {
  if(.Toeplitz_hasAcf(object$cpp_ptr)) {
    obj.acf <- object$getAcf()[1:min(6, object$size)]
    obj.acf <- signif(obj.acf, digits = 3)
    if(object$size > 6) obj.acf <- c(obj.acf, "...")
  } else {
    obj.acf <- "NULL"
  }
  cat("Toeplitz matrix of size", object$size, "\n",
      "acf: ", obj.acf, "\n")
})

# ncol
#' @export
setMethod("ncol", "Toeplitz", function(x){
  x$size
})

# nrow
#' @export
setMethod("nrow", "Toeplitz", function(x){
  x$size
})

# dim
#' @export
setMethod("dim", "Toeplitz", function(x){
  rep(x$size, 2)
})

# Toeplitz-Matrix multiplication
#' @export
setMethod("%*%", signature(x = "Toeplitz", y = "ANY"), function(x, y) {
  if(!.Toeplitz_hasAcf(x$cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(is.vector(y)) y <- as.matrix(y)
  if(!is.matrix(y)) {
    stop("Second argument should be a matrix or vector.")
  }
  if(nrow(y) != x$size) {
    stop("Toeplitz and second argument are non-conformable.")
  }
  .Toeplitz_Multiply(x$cpp_ptr, y)
})

# Matrix-Toeplitz multiplication
setMethod("%*%", signature(x = "ANY", y = "Toeplitz"), function(x, y) {
  if(!.Toeplitz_hasAcf(y$cpp_ptr)) {
    stop("setAcf has not been called yet")
  }
  if(is.vector(x)) x <- as.matrix(x)
  if(!is.matrix(x)) {
    stop("First argument should be a matrix or vector.")
  }
  if(ncol(x) != y$size) {
    stop("First argument and Toeplitz are non-conformable.")
  }
  t(.Toeplitz_Multiply(y$cpp_ptr, t(x)))
})

# determinant
#' @export
setMethod("determinant", "Toeplitz",
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
setMethod("solve", "Toeplitz", function(a, b, ...) {
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

# now make methods display arguments
.DollarNames.Toeplitz <- function(x, pattern)
    grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
