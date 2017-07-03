#' @title Convert position to increment autocorrelations.
#'
#' @description Converts the autocorrelation of a stationary sequence \code{X} to that of its increments, \code{dX == diff(X)}.
#' @param acf Length-\code{N} vector of position autocorrelations.
#' @return Length \code{N-1} vector if increment autocorrelations.
#' @examples
#' acf2incr(acf = exp(-(0:10)))
#' @export
acf2incr <- function(acf) {
  N <- length(acf)-1
  if(N == 1) {
    iacf <- 2*(acf[1]-acf[2])
  } else {
    iacf <- 2*acf[1:N] - acf[1:N+1] - acf[c(2, 1:(N-1))]
  }
  iacf
}
