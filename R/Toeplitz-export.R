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
Toep.acf <- function(Toeplitz, acf){}

#' @export
traceT2 <- function(Toeplitz, acf1){}

#' @export
traceT4 <- function(Toeplitz, acf1, acf2){}