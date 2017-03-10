loadModule("Class_Toeplitz", TRUE)

setClass("Toeplitz_Cpp")
setRcppClass(Class = "Toeplitz_Cpp",
             CppClass = "Toeplitz_Cpp",
             module = "Class_Toeplitz",
             saveAs = ".Toeplitz")

setMethod("Toep.det", "Toeplitz_Cpp", function(Toeplitz, logarithm = TRUE) {
  ldT <- Toeplitz$Det()
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


setMethod("Toep.solve", "Toeplitz_Cpp", function(Toeplitz, x){
  if(!(is.matrix(x) || is.vector(x))){
    stop("argument x should be either matrix or vector")
  }

  if(is.matrix(x)){
    if(nrow(x) == Toeplitz$DimCheck()){
      mat <- Toeplitz$Solve(x)
    }
    else{
      stop("incompatible dimension of x")
    }
  }
  if(is.vector(x)){
    if(length(x) == Toeplitz$DimCheck()){
      mat <- Toeplitz$SolveVec(x)
    }
    else{
      stop("incompatible dimension of x")
    }
  }
  mat
})


setMethod("Toep.trace", "Toeplitz_Cpp", function(Toeplitz, acf1){
  if(!(is.matrix(acf1) || is.vector(acf1))){
    stop("argument acf1 should be either matrix or vector")
  }
  if(length(acf1) != Toeplitz$DimCheck()){
    stop("incompatible dimension of acf1")
  }
  Toeplitz$TraceProd(acf1)
})

setMethod("Toep.deriv", "Toeplitz_Cpp", function(Toeplitz, acf1, acf2){
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