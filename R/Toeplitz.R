#' @export determinant
#' @export solve
#' @export nrow
#' @export ncol
NULL

#' Creating a Toeplitz class
setClass("Toeplitz_Cpp")
setRcppClass(Class = "Toeplitz_Cpp",
             CppClass = "Toeplitz_Cpp",
             module = "Class_Toeplitz",
             saveAs = ".Toeplitz")

#' Creating the Determinant function
setMethod("determinant", "Toeplitz_Cpp", function(x, logarithm) {
  ldT <- x$Det()
  if(!logarithm){
    ldT <- exp(ldT)
  }
  ldT
})

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

#' display method
setMethod("show", "Toeplitz_Cpp", function(object){
  if(object$flag_acf()){
    obj.acf <- c(round(object$acf[1:6], digits = 3), "...")
  }else{
    obj.acf <- "NULL"
  }
  cat("Toeplitz object of size", object$DimCheck(), "\n", 
      "First row of Toeplitz matrix: ", obj.acf)
})

#' number of col of Toeplitz matrix
setMethod("ncol", "Toeplitz_Cpp", function(x){
  x$DimCheck()
})

#' number of row of Toeplitz matrix
setMethod("nrow", "Toeplitz_Cpp", function(x){
  x$DimCheck()
})

#' dimension of Toeplitz matrix
setMethod("dim", "Toeplitz_Cpp", function(x){
  x$DimCheck()
})

#' Toeplitz matrix times vector(matrix)
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

#' Inversion of Toeplitz matrix times vector(matrix)
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

#' # #' @export
#' # Toep.acf <- function(Toeplitz, acf){}
#' 
#' #' Creating the interface for inputing acf
#' setGeneric("Toep.acf", function(Toeplitz, acf){})
#' setMethod("Toep.acf", "Toeplitz_Cpp", function(Toeplitz, acf){
#'   if(!(is.vector(acf) || is.matrix(acf))){
#'     stop("argument acf should be either matrix or vector")
#'   }
#'   if(length(acf) != Toeplitz$DimCheck()){
#'     stop("incompatible dimension of acf")
#'   }
#'   Toeplitz$AcfInput(acf)
#' })

#' # #' @export
#' # traceT2 <- function(Toeplitz, acf1){}
#' 
#' #' trace of solve(T1, T2), as T1, T2 Toeplitz
#' setGeneric("traceT2", function(Toeplitz, acf1){})
#' setMethod("traceT2", "Toeplitz_Cpp", function(Toeplitz, acf1){
#'   if(!(is.matrix(acf1) || is.vector(acf1))){
#'     stop("argument acf1 should be either matrix or vector")
#'   }
#'   if(length(acf1) != Toeplitz$DimCheck()){
#'     stop("incompatible dimension of acf1")
#'   }
#'   Toeplitz$TraceProd(acf1)
#' })
#' 
#' # #' @export
#' # traceT4 <- function(Toeplitz, acf1, acf2){}
#' 
#' #' trace of solve(T1, T2) %*% solve(T1, T3), as T1 T2 T3 Toeplitz
#' setGeneric("traceT4", function(Toeplitz, acf1, acf2){})
#' setMethod("traceT4", "Toeplitz_Cpp", function(Toeplitz, acf1, acf2){
#'   if(!(is.matrix(acf1) || is.vector(acf1))){
#'     stop("argument acf1 should be either matrix or vector")
#'   }
#'   if(!(is.matrix(acf2) || is.vector(acf2))){
#'     stop("argument acf2 should be either matrix or vector")
#'   }
#'   if(length(acf1) != Toeplitz$DimCheck()){
#'     stop("incompatible dimension of x")
#'   }
#'   if(length(acf2) != Toeplitz$DimCheck()){
#'     stop("incompatible dimension of y")
#'   }
#'   Toeplitz$TraceDeriv(acf1, acf2)
#' })