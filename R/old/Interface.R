## Interface for SuperGauss
## Generic methods for RC class

require(Rcpp)
require(SuperGauss)

n <- 10
acf <- exp(-(1:10))
X <- 1:10
Y <- matrix(1, 10, 1)

## # class definition
## setRcppClass(Class = ".Toeplitz", CppClass = "Toeplitz",
##              module = "Toeplitz", saveAs = ".Toeplitz")

## # constructor
## ToeplitzCons <- function(n, acf) {
##   if(missing(acf)) {
##     if(missing(n)) {
##       stop("Must provide n or acf.")
##     } else {
##       if(length(n) != 1) stop("n must be an integer.")
##       T <- new(Toeplitz, n)
##     }
##   } else {
##     n <- length(acf)
##     T <- new(Toeplitz, n)
##     T$AcfInput(acf)
##   }
##   T
## }

T1 <- Toeplitz(n = length(acf))

solve(T1, Y)
## setMethod("solve", "Toeplitz_CPP", function(a, b) {
##   a$Solve(b)
## })

Toeplitz_CPP <- setRcppClass(Class = "Toeplitz_CPP")

.Toeplitz <- setRefClass(Class = "Toeplitz",
                         contains = "Rcpp_Toeplitz_CPP",
                         methods = list(initialize = function(n, acf) {
                           Toeplitz(n, acf)
                         }))

setMethod("solve", "Toeplitz", function(a, b) {
  if(missing(b)) b <- diag(a$DimCheck())
  a$Solve(b)
})


# generic show method
setMethod("show", "Toeplitz_CPP", function(object){
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
