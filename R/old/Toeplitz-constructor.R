#' Toeplitz Class Constructor.
#'
#' @export
ToeplitzConst <- function(n, acf) {
  if(missing(acf)) {
    if(missing(n)) {
      stop("Must provide n or acf.")
    } else {
      if(length(n) != 1) stop("n must be an integer.")
      ## T <- new(Toeplitz_CPP, n)
      T <- .Toeplitz(n)
    }
  } else {
    n <- length(acf)
    T <- new(Toeplitz_CPP, n)
    T$AcfInput(acf)
  }
  T
}

#' Toeplitz Class
#'
#' @export Toeplitz
#' @exportClass Toeplitz
setRcppClass(Class = "Toeplitz", CppClass = "Toeplitz_CPP",
             module = "Toeplitz",
             saveAs = ".Toeplitz")

setMethod("solve", "Toeplitz", function(a, b) {
  if(missing(b)) b <- diag(a$DimCheck())
  a$Solve(b)
})
