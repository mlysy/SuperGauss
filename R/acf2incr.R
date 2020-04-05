#' Convert position autocorrelations to increment autocorrelations.
#'
#' Convert the autocorrelation of a stationary sequence `x = (x_1, ..., x_N)` to that of its increments, `dx = (x_2 - x_1, ..., x_N - x_(N-1))`.
#'
#' @param acf Length-`N` vector of position autocorrelations.
#' @return Length `N-1` vector of increment autocorrelations.
#'
#' @example examples/acf2incr.R
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
