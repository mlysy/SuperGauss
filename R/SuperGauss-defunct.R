#' Defunct functions in \pkg{SuperGauss}.
#'
#' @name SuperGauss-defunct
#' @section The following functions have been removed from the \pkg{SuperGauss} package:
#' \describe{
#'   \item{`rSnorm()`}{Please use [rnormtz()] instead.}
#'   \item{`dSnorm()`}{Please use [dnormtz()] instead.}
#'   \item{`Snorm.grad()`}{Please use the `grad()` method in the [NormalToeplitz] class.}
#'   \item{`Snorm.hess()`}{Please use the `hess()` method in the [NormalToeplitz] class.}
#' }
NULL

#' @rdname SuperGauss-defunct
#' @usage NULL
#' @export
rSnorm <- function(...) {
  .Defunct(new = "rnormtz", package= "SuperGauss")
}

#' @rdname SuperGauss-defunct
#' @usage NULL
#' @export
dSnorm <- function(...) {
  .Defunct(new = "dnormtz", package= "SuperGauss")
}

#' @rdname SuperGauss-defunct
#' @usage NULL
#' @export
Snorm.grad <- function(...) {
  .Defunct(new = "NormalToeplitz$grad", package= "SuperGauss")
}

#' @rdname SuperGauss-defunct
#' @usage NULL
#' @export
Snorm.hess <- function(...) {
  .Defunct(new = "NormalToeplitz$hess", package= "SuperGauss")
}
