#' Radial basis kernel function
#' @param lag Numeric. Lags at which the kernel function is evaluated.
#' @param var Scalar widthscale parameter.
#' @param ell Scalar lengthscale parameter.
#' @return Numeric of the same length as `lag`.
#' @export
rbf_acf <- function(lag, var=1, ell=0.1){
  var * exp(- abs(lag)^2 / ell)
}

#' Long memory latent process ARFIMA(0, d, 0)
arfima_acf <- function(lag, var, d){
  delta_t <- lag[2]-lag[1]
  n <- lag + delta_t
  var*(n)^{2*d-1} 
}

#' Long memory rational quadratic kernel rquad(alpha, k)
#' @param k Lengthscale.
#' @param d Additional tuning parameter which has a 1-to-1 transform with alpha.
#' @param alpha Only one of alpha and d should be provided. 
#' @details
#' When k=1 and as lag goes to infinity, its acf has the same asymptotic 
#' behavior as ARIMA(0, d, 0). Note that d should satisfy |d| < 0.5, see 
#' Robinson (1995) Log-Periodogram Regression of Time Series with Long Range 
#' Dependence, AOS.
rquad_acf <- function(lag, var, k, d, alpha=NULL){
  if (is.null(alpha)){
    alpha = 0.5-d
  } else if (is.numeric(d)){
      stop("Only one of d and alpha should be provided.")
  }
  delta_t <- lag[2]-lag[1]
  var * (1 + lag^2/(alpha*k^2))^(-alpha)
}
