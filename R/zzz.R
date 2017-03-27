loadModule("Class_Toeplitz", TRUE)

setClass("Toeplitz_Cpp")
setRcppClass(Class = "Toeplitz_Cpp",
             CppClass = "Toeplitz_Cpp",
             module = "Class_Toeplitz",
             saveAs = ".Toeplitz")

setMethod("determinant", "Toeplitz_Cpp", function(x, logarithm, ...) {
  ldT <- x$Det()
  if(!logarithm){
    ldT <- exp(ldT)
  }
  ldT
})

setMethod("Toep.acf", "Toeplitz_Cpp", function(Toeplitz, acf){
  if(!(is.vector(acf) || is.matrix(acf))){
    stop("argument acf should be either matrix or vector")
  }
  if(length(acf) != Toeplitz$DimCheck()){
    stop("incompatible dimension of acf")
  }
  Toeplitz$AcfInput(acf)
})


setMethod("show", "Toeplitz_Cpp", function(object){
  print(object$AcfCheck())
})


setMethod("dim", "Toeplitz_Cpp", function(x){
  x$DimCheck()
})

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


setMethod("solve", "Toeplitz_Cpp", function(a, b, ...){
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


setMethod("traceT2", "Toeplitz_Cpp", function(Toeplitz, acf1){
  if(!(is.matrix(acf1) || is.vector(acf1))){
    stop("argument acf1 should be either matrix or vector")
  }
  if(length(acf1) != Toeplitz$DimCheck()){
    stop("incompatible dimension of acf1")
  }
  Toeplitz$TraceProd(acf1)
})

setMethod("traceT4", "Toeplitz_Cpp", function(Toeplitz, acf1, acf2){
  if(!(is.matrix(acf1) || is.vector(acf1))){
    stop("argument acf1 should be either matrix or vector")
  }
  if(!(is.matrix(acf2) || is.vector(acf2))){
    stop("argument acf2 should be either matrix or vector")
  }
  if(length(acf1) != Toeplitz$DimCheck()){
    stop("incompatible dimension of x")
  }
  if(length(acf2) != Toeplitz$DimCheck()){
    stop("incompatible dimension of y")
  }
  Toeplitz$TraceDeriv(acf1, acf2)
})