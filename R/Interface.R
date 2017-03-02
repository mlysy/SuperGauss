## Interface for SuperGauss
## Generic methods for RC class

require(SuperGauss)

n <- 10
acf <- 1:10
X <- 1:10
Y <- matrix(1, 10, 1)

# class definition
.Toeplitz <- setRefClass("ToeplitzMat", fields = list(toep = "Rcpp_Toeplitz"))

# constructor
ToeplitzMat <- function(acf){
  n <- length(acf)
  toep <- new(Toeplitz, n)
  toep$AcfInput(acf)
  .Toeplitz$new(toep = toep)
}

T1 <- ToeplitzMat(acf)

# generic show method
setMethod("show", "ToeplitzMat", function(object){
  .toep <- object$toep
  print(.toep$DimCheck())
})

T1

# generic product method
setMethod("%*%", c(x = "ToeplitzMat", y = "matrix"), function(x, y){
  .toep <- x$toep
  vec <- .toep$Mult(y)
  vec
})

T1 %*% Y

# generic solve method
setMethod("solve", "ToeplitzMat", function(a, b){
  .toep <- a$toep
  .toep$Solve(b)
})

# generic log-determinant method
setMethod("Toeplitz.Det", "ToeplitzMat", function(obj){
  .toep <- obj$toep
  .toep$Det()
})

# generic trace-prod method
setMethod("Toeplitz.Prod.Trace", "ToeplitzMat", function(obj, acf2){
  .toep <- obj$toep
  .toep$TraceProd(acf2)
})

# generic trace-deriv method
setMethod("Toeplitz.Deriv.Trace", "ToeplitzMat", function(obj, acf2, acf3){
  .toep <- obj$toep
  .toep$TraceDriv(acf2, acf3)
})