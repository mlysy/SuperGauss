#' Toeplitz Matrix Class
#'
#' Contains methods for efficient linear algebra involving symmetric positive definite Toeplitz matrices.
#' @details Toeplitz objects are created with \code{Toeplitz}
#' @name Toeplitz-Class
#' @export determinant
#' @export solve
#' @export nrow
#' @export ncol
NULL

# Creating a Toeplitz class
setClass("Toeplitz_Cpp")
setRcppClass(Class = "Toeplitz_Cpp",
             CppClass = "Toeplitz_Cpp",
             module = "Class_Toeplitz",
             saveAs = ".Toeplitz")




#' @rdname Toeplitz-Class
#' @usage x %*% y
#' @param y A regular matrix or vector.
setMethod("%*%", "Toeplitz_Cpp", function(x, y){
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

# display method
setMethod("show", "Toeplitz_Cpp", function(object){
  if(object$flag_acf()){
    obj.acf <- c(round(object$acf[1:6], digits = 3), "...")
  }else{
    obj.acf <- "NULL"
  }
  cat("Toeplitz object of size", object$DimCheck(), "\n",
      "First row of Toeplitz matrix: ", obj.acf)
})



#' @rdname Toeplitz-Class
#' @usage determinant(T, logarithm = TRUE, ...)
#' @param x A \code{Toeplitz_Matrix} object.
#' @param logarithm Logical; whether the determinant should be calculated on the log scale.
setMethod("determinant", "Toeplitz_Cpp",
          function(x, logarithm = TRUE, ...) {
  ldT <- x$Det()
  if(!logarithm){
    ldT <- exp(ldT)
  }
  ldT
})


#' @rdname Toeplitz-Class
#' @usage solve(a, b)
#' @param a A \code{Toeplit_Matrix} object.
#' @param b A regular matrix or vector.
setMethod("solve", "Toeplitz_Cpp", function(a, b){
  if(!(is.matrix(b) || is.vector(b))){
    stop("argument x should be either matrix or vector")
  }
  if(is.matrix(b)){
    if(nrow(b) == a$DimCheck()){
      mat <- a$Solve(b)
    }
    else{
      stop("incompatible dimension of b")
    }
  }
  if(is.vector(b)){
    if(length(b) == a$DimCheck()){
      mat <- a$SolveVec(b)
    }
    else{
      stop("incompatible dimension of b")
    }
  }
  mat
})

# number of col of Toeplitz matrix
setMethod("ncol", "Toeplitz_Cpp", function(x){
  x$DimCheck()
})

# number of row of Toeplitz matrix
setMethod("nrow", "Toeplitz_Cpp", function(x){
  x$DimCheck()
})

# dimension of Toeplitz matrix
setMethod("dim", "Toeplitz_Cpp", function(x){
  rep(x$DimCheck(), 2)
})

#' @rdname Toeplitz-Class
#' @param n Dimension of Toeplitz matrix.
#' @param acf Autocorrelation of Toeplitz matrix, i.e., first row or column.
#' @export
Toeplitz <- function(n, acf){
  if(missing(n)){
    n <- length(acf)
  }
  T <- .Toeplitz$new(n)
  if(missing(acf)){
    message("please use Toeplitz$setAcf to input the acf of Toeplitz variance")
  } else{
    T$AcfInput(acf)
  }
  T
}
