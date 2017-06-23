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

#' Toeplitz Matrix Class
#'
#' Contains methods for efficient linear algebra involving Toeplitz variance matrices.
#' @details \code{Toeplitz_Matrix} objects are created with the R constructor \code{Toeplitz}, which accepts an autocorrelation vector \code{acf}, or an integer \code{n} which can be used to preallocate memory.  The \code{acf} of the \code{Toeplitz_Matrix} is accessed with \code{setAcf} and \code{getAcf}.
#' The \code{Toeplitz_Matrix} class efficiently implements some generic methods Toeplitz variance matrix calculations, e.g., %*% (multiplication), \code{solve}, \code{determinant}.  Non-generic methods implemented for Toeplitz matrices are:
#' \itemize{
#'   \item \code{setAcf(acf1)}: input the acf of Toeplitz variance
#'   \item \code{getAcf()}: obtain the inputed acf
#'   \item \code{traceT2(acf2)}: matrix trace of \code{solve(.self, toeplitz(acf2))}.
#'   \item \code{traceT4(acf2, acf3)}: matrix trace of \code{solve(.self, toeplitz(acf2))} %*% \code{solve(.self, toeplitz(acf3))}.
#' }
#' These are used for calculating the gradient and Hessian of multivariate normal loglikelihoods having Toeplitz variance matrix.
#' @name Toeplitz-Class
#' @export determinant
#' @export solve
#' @export nrow
#' @export ncol
NULL


# class skeleton
.Toeplitz <- setRefClass("Toeplitz_Matrix",
                         fields = list(cpp_ptr = "externalptr",
                                       size = "numeric"))
.Toeplitz$lock("cpp_ptr")
.Toeplitz$lock("size")

# non-generic methods

# internal constructor
.Toeplitz$methods(initialize = function(n) {
  cpp_ptr <<- .Toeplitz_constructor(n)
  size <<- n
})
# getter/setter
.Toeplitz$methods(setAcf = function(acf) {
  if(length(acf) != size) {
    stop("acf has wrong length.")
  }
  .Toeplitz_setAcf(cpp_ptr, acf)
})
.Toeplitz$methods(getAcf = function() {
  .Toeplitz_getAcf(cpp_ptr)
})
# traceT2/traceT4
.Toeplitz$methods(traceT2 = function(acf2) {
  if(length(acf2) != size) {
    stop("acf2 has wrong length.")
  }
  .Toeplitz_traceT2(cpp_ptr, acf2)
})
.Toeplitz$methods(traceT4 = function(acf2, acf3) {
  if(length(acf2) != size) {
    stop("acf2 has wrong length.")
  }
  if(length(acf3) != size) {
    stop("acf3 has wrong length.")
  }
  .Toeplitz_traceT4(cpp_ptr, acf2, acf3)
})

# generic methods

# show method
setMethod("show", "Toeplitz_Matrix", function(object) {
  if(.Toeplitz_hasAcf(object$cpp_ptr)) {
    obj.acf <- object$getAcf()[1:min(6, object$size)]
    obj.acf <- round(obj.acf, digits = 3)
    if(object$size > 6) obj.acf <- c(obj.acf, "...")
  } else {
    obj.acf <- "NULL"
  }
  cat("Toeplitz_Matrix of size", object$size, "\n",
      "acf: ", obj.acf, "\n")
})

#' @rdname Toeplitz-Class
#' @usage ncol(x)
#' @examples
#' ncol(Toep)
setMethod("ncol", "Toeplitz_Matrix", function(x){
  x$size
})

#' @rdname Toeplitz-Class
#' @usage nrow(x)
#' @examples
#' nrow(Toep)
setMethod("nrow", "Toeplitz_Matrix", function(x){
  x$size
})

#' @rdname Toeplitz-Class
#' @usage dim(x)
#' @examples
#' dim(Toep)
setMethod("dim", "Toeplitz_Matrix", function(x){
  rep(x$size, 2)
})


#' @rdname Toeplitz-Class
#' @usage x %*% y
#' @param y A regular matrix or vector.
#' @examples
#' Y <- matrix(rnorm(N*2), N, 2)
#' Toep %*% Y
setMethod("%*%", signature(x = "Toeplitz_Matrix", y = "ANY"), function(x, y) {
  if(is.vector(y)) y <- as.matrix(y)
  if(!is.matrix(y)) {
    stop("Second argument should be a matrix or vector.")
  }
  if(nrow(y) != x$size) {
    stop("Toeplitz_Matrix and second argument are non-conformable.")
  }
  .Toeplitz_Multiply(x$cpp_ptr, y)
})
setMethod("%*%", signature(x = "ANY", y = "Toeplitz_Matrix"), function(x, y) {
  if(is.vector(x)) x <- as.matrix(x)
  if(!is.matrix(x)) {
    stop("First argument should be a matrix or vector.")
  }
  if(nrow(y) != x$size) {
    stop("First argument and Toeplitz_Matrix are non-conformable.")
  }
  t(.Toeplitz_Multiply(y$cpp_ptr, t(x)))
})

#' @rdname Toeplitz-Class
#' @usage determinant(x, logarithm = TRUE, ...)
#' @param x A \code{Toeplitz_Matrix} object.
#' @param logarithm Logical; whether the determinant should be calculated on the log scale.
#' @examples
#' determinant(Toep, logarithm = TRUE)
setMethod("determinant", "Toeplitz_Matrix",
          function(x, logarithm = TRUE, ...) {
  ldT <- .Toeplitz_Determinant(x$cpp_ptr)
  if(!logarithm) {
    ldT <- exp(ldT)
  }
  ldT
})

#' @rdname Toeplitz-Class
#' @usage solve(a, b)
#' @param a A \code{Toeplitz_Matrix} object.
#' @param b A regular matrix or vector.
#' @examples
#' Y <- matrix(rnorm(N*2), N, 2)
#' solve(Toep, Y)
setMethod("solve", "Toeplitz_Matrix", function(a, b){
  if(is.vector(b)) b <- as.matrix(b)
  if(!is.matrix(b)) {
    stop("b must be a matrix or vector.")
  }
  if(nrow(b) != a$size) {
    stop("a and b are non-conformable.")
  }
  .Toeplitz_Solve(a$cpp_ptr, b)
})

#' @rdname Toeplitz-Class
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
