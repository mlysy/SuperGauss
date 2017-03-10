#' @export
Toeplitz <- function(n, acf){
  if(missing(n)){
    n <- length(acf)
  }
  T <- .Toeplitz$new(n)
  if(missing(acf)){
    message("please use $AcfInput to set the acf of Toeplitz function")
  } else{
    T$AcfInput(acf)
  }
  T
}

#' @export
Toep.det <- function(Toeplitz, logarithm = TRUE){}

#' @export
show <- function(object){}

#' @export
Toep.acf <- function(Toeplitz, acf){}

#' @export
Toep.solve <- function(Toeplitz, x){}

#' @export
dim <- function(x){}

#' @export
Toep.trace <- function(Toeplitz, acf1){}

#' @export
Toep.deriv <- function(Toeplitz, acf1, acf2){}