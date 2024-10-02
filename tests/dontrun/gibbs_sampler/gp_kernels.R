#' Radial basis kernel function
#' @param lag Numeric. Lags at which the kernel function is evaluated.
#' @param var Scalar widthscale parameter.
#' @param ell Scalar lengthscale parameter.
#' @return Numeric of the same length as `lag`.
#' @export
rbf_acf <- function(lag, var=1, ell=0.1){
  var * exp(- abs(lag)^2 / ell)
}

# Long memory latent process ARFIMA(0, d, 0)
arfima_acf <- function(lag, var, d){
  delta_t <- lag[2]-lag[1]
  n <- lag + delta_t
  var*(n)^{2*d-1} 
}
# Long memory rational quadratic kernel rquad(alpha, k)
rquad_acf <- function(lag, var, k, alpha){
  delta_t <- lag[2]-lag[1]
  var * (1 + lag^2/(alpha*k^2))^(-alpha)
}