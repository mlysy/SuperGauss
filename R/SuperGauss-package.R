#' Superfast inference for stationary Gaussian time series.
#'
#' @details While likelihood calculations with stationary Gaussian time series generally scale as `O(N^2)` in the number of observations, this package implements an algorithm which scales as `O(N log^2 N)`.  "Superfast" algorithms for loglikelihood gradients and Hessians are also provided.  The underlying C++ code is distributed through a header-only library found in the installed package's `include` directory.
#'
#' @example examples/SuperGauss-package.R
#' @importFrom Rcpp evalCpp
#' @importFrom methods new show
#' @importFrom stats rnorm
#' @import R6
#' @useDynLib SuperGauss, .registration = TRUE
"_PACKAGE"
