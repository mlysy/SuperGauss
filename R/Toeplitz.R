#' Toeplitz Matrix Class
#'
#' Contains methods for efficient linear algebra involving Toeplitz variance matrices.
#' @details \code{Toeplitz_Matrix} objects are created with the R constructor \code{Toeplitz}, which accepts an autocorrelation vector \code{acf}, or an integer \code{n} which can be used to preallocate memory.  The \code{acf} of the \code{Toeplitz_Matrix} is accessed with \code{setAcf} and \code{getAcf}.
#' The \code{Toeplitz_Matrix} class efficiently implements some generic methods Toeplitz variance matrix calculations, e.g., %*% (multiplication), \code{solve}, \code{determinant}.  Non-generic methods implemented for Toeplitz matrices are:
#' \itemize{
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

# Creating a Toeplitz class
setClass("Toeplitz_Matrix") # otherwise install has benign error message
setRcppClass(Class = "Toeplitz_Matrix",
             CppClass = "Toeplitz_Cpp",
             module = "Toeplitz_Class",
             saveAs = ".Toeplitz")


#' @rdname Toeplitz-Class
#' @usage .self$setAcf(acf)
#' @usage .self$getAcf()
#' @param n Integer dimension of Toeplitz variance matrix.
#' @param acf Autocorrelation vector of Toeplitz matrix, i.e., first row or column.
#' @export
Toeplitz <- function(n, acf){
  if(missing(n)){
    n <- length(acf)
  }
  Tz <- .Toeplitz$new(n)
  if(!missing(acf)) {
    Tz$setAcf(acf)
  }
  Tz
}

# show method
setMethod("show", "Toeplitz_Matrix", function(object){
  if(object$flag_acf()){
    obj.acf <- object$getAcf()[1:min(6, object$DimCheck())]
    obj.acf <- round(obj.acf, digits = 3)
    if(length(obj.acf) > 6) obj.acf <- c(obj.acf, "...")
  }else{
    obj.acf <- "NULL"
  }
  cat("Toeplitz_Matrix of size", object$DimCheck(), "\n",
      "ACF: ", obj.acf, "\n")
})


#' @rdname Toeplitz-Class
#' @usage ncol(x)
setMethod("ncol", "Toeplitz_Matrix", function(x){
  x$DimCheck()
})

#' @rdname Toeplitz-Class
#' @usage nrow(x)
setMethod("nrow", "Toeplitz_Matrix", function(x){
  x$DimCheck()
})

#' @rdname Toeplitz-Class
#' @usage dim(x)
setMethod("dim", "Toeplitz_Matrix", function(x){
  rep(x$DimCheck(), 2)
})


#' @rdname Toeplitz-Class
#' @usage x %*% y
#' @param y A regular matrix or vector.
setMethod("%*%", "Toeplitz_Matrix", function(x, y){
  if(!(is.matrix(y) || is.vector(y))){
    stop("argument y should be either matrix or vector")
  }
  if(is.matrix(y)){
    if(nrow(y) == x$DimCheck()){
      mat <- x$Mult(y)
    }
    else{
      stop("incompatible dimension of y")
    }
  }
  if(is.vector(y)){
    if(length(y) == x$DimCheck()){
      mat <- x$MultVec(y)
    }
    else{
      stop("incompatible dimension of y")
    }
  }
  mat
})

#' @rdname Toeplitz-Class
#' @usage determinant(x, logarithm = TRUE, ...)
#' @param x A \code{Toeplitz_Matrix} object.
#' @param logarithm Logical; whether the determinant should be calculated on the log scale.
setMethod("determinant", "Toeplitz_Matrix",
          function(x, logarithm = TRUE, ...) {
  ldT <- x$Det()
  if(!logarithm){
    ldT <- exp(ldT)
  }
  ldT
})


#' @rdname Toeplitz-Class
#' @usage solve(a, b)
#' @param a A \code{Toeplitz_Matrix} object.
#' @param b A regular matrix or vector.
setMethod("solve", "Toeplitz_Matrix", function(a, b){
  if(!(is.matrix(b) || is.vector(b))){
    stop("b must be a matrix or vector.")
  }
  if(is.matrix(b)){
    if(nrow(b) == a$DimCheck()){
      mat <- a$Solve(b)
    }
    else{
      stop("a and b have incompatible dimensions.")
    }
  }
  if(is.vector(b)){
    if(length(b) == a$DimCheck()){
      mat <- a$SolveVec(b)
    }
    else{
      stop("a and b have incompatible dimensions.")
    }
  }
  mat
})
