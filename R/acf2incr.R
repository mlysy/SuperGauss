#' @title Position autocorrelation(ACF) to Increment ACF.
#'
#' @description Converts the ACF of a stationary sequence \code{{X0, X1, ..., XN}} to those of
#' its increments \code{{dX1, dX2, ..., dXN}}.
#' @param gam length \code{N} vector of position ACF
#' @return Length \code{N-1} vector if increment ACF.
#' @examples
#' acf2incr(rnorm(10))
#' @export
acf2incr <- function(gam) {
  N <- length(gam)-1
  if(N == 1) {
    igam <- 2*(gam[1]-gam[2])
  } else {
    igam <- 2*gam[1:N] - gam[1:N+1] - gam[c(2, 1:(N-1))]
  }
  igam
}