#' Multivariate normal with Toeplitz variance matrix.
#'
#' @description Provides methods for the Normal-Toeplitz (NTz) distribution defined as
#' ```
#' z ~ NTz(acf)   <=>   z ~ Normal(0, toeplitz(acf)),
#' ```
#' i.e., for a multivariate normal with mean zero and variance `Tz = toeplitz(acf)`.
#' @export
NormalToeplitz <- R6Class(
  classname = "NormalToeplitz",

  private = list(

    NTz_ptr = NULL,
    size = NA

  ),

  public = list(

    #' @description Class constructor.
    #'
    #' @param n Size of the Toeplitz matrix.
    #' @return A `NormalToeplitz` object.
    initialize = function(n) {
      private$NTz_ptr <- .NormalToeplitz_constructor(n)
      private$size <- n
    },

    #' @description Log-density function.
    #'
    #' @param z Density argument. A vector of length `n` or an `n x N` matrix where each column is an observation.
    #' @param acf A vector of length `n` containing the first row/column of the Toeplitz variance matrix.
    #' @return A scalar or vector of length `N` containing the log-density of the NTz evaluated at its arguments.
    logdens = function(z, acf) {
      check_var(z, N = private$size, varname = "z")
      check_var(acf, N = private$size, varname = "acf")
      .NormalToeplitz_logdens(private$NTz_ptr, z, acf)
    },

    #' @description Gradient of the log-density with respect to parameters.
    #'
    #' @param z Density argument.  A vector of length `n`.
    #' @param dz An `n x ntheta` matrix containing the gradient `dz/dtheta`.
    #' @param acf A vector of length `n` containing the first row/column of the Toeplitz variance matrix.
    #' @param dacf An `n x ntheta` matrix containing the gradient `dacf/dtheta`.
    #' @return A vector of length `ntheta` containing the gradient of the NTz log-density with respect to `theta`.
    grad = function(z, dz, acf, dacf) {
      check_var(z, dx = dz, N = private$size, varname = "z")
      check_var(acf, dx = dacf, N = private$size, varname = "acf")
      if(ncol(dz) != ncol(dacf)) {
        stop("dz and dacf must have the same number of columns.")
      }
      ntheta <- ncol(dz)
      .NormalToeplitz_grad(private$NTz_ptr, z, dz, acf, dacf, ntheta)
    },

    #' @description Hessian of log-density with respect to parameters.
    #'
    #' @param z Density argument.  A vector of length `n`.
    #' @param dz An `n x ntheta` matrix containing the gradient `dz/dtheta`.
    #' @param d2z An `n x ntheta x ntheta` array containing the Hessian `d^2z/dtheta^2`.
    #' @param acf A vector of length `n` containing the first row/column of the Toeplitz variance matrix.
    #' @param dacf An `n x ntheta` matrix containing the gradient `dacf/dtheta`.
    #' @param d2acf An `n x ntheta x ntheta` array containing the Hessian `dacf^2/dtheta^2`.
    #' @return An `ntheta x ntheta` matrix containing the Hessian of the NTz log-density with respect to `theta`.
    hess = function(z, dz, d2z, acf, dacf, d2acf) {
      check_var(z, dx = dz, d2x = d2z, N = private$size, varname = "z")
      check_var(acf, dx = dacf, d2x = d2acf, N = private$size,
                varname = "acf")
      if(ncol(dz) != ncol(dacf)) {
        stop("dz and dacf must have the same number of columns.")
      }
      if(!all(dim(d2z) == dim(d2acf))) {
        stop("d2z and d2acf must have the same dimensions.")
      }
      ntheta <- ncol(dz)
      .NormalToeplitz_hess(private$NTz_ptr, z, dz,
                           matrix(d2z, private$size, ntheta*ntheta),
                           acf, dacf,
                           matrix(d2acf, private$size, ntheta*ntheta),
                           ntheta)
    },

    #' @description Full gradient of log-density function.
    #'
    #' @param z Density argument.  A vector of length `n`.
    #' @param acf A vector of length `n` containing the first row/column of the Toeplitz variance matrix.
    #' @return A `n x 2` matrix containing the gradient of the NTz log-density with respect to each element of `z` and `acf`.
    grad_full = function(z, acf) {
      check_var(z, N = private$size, varname = "z")
      check_var(acf, N = private$size, varname = "acf")
      .NormalToeplitz_grad_full(private$NTz_ptr, z, acf)
    }
  )

)

#--- custom methods ------------------------------------------------------------

## # Class constructor
## NormalToeplitz$set("public", "initialize", overwrite = TRUE, function(n) {
##   private$NTz_ptr <- .NormalToeplitz_constructor(n)
##   private$size <- n
## })

## # Log-density function.
## NormalToeplitz$set("public", "logdens", overwrite = TRUE, function(z, acf) {
##   check_var(z, N = private$size, varname = "z")
##   check_var(acf, N = private$size, varname = "acf")
##   .NormalToeplitz_logdens(private$NTz_ptr, z, acf)
## })

## # Gradient of log-density with respect to parameters.
## NormalToeplitz$set("public", "grad", overwrite = TRUE, function(z, dz, acf, dacf) {
##   check_var(z, dx = dz, N = private$size, varname = "z")
##   check_var(acf, dx = dacf, N = private$size, varname = "acf")
##   if(ncol(dz) != ncol(dacf)) {
##     stop("dz and dacf must have the same number of columns.")
##   }
##   ntheta <- ncol(dz)
##   .NormalToeplitz_grad(private$NTz_ptr, z, dz, acf, dacf, ntheta)
## })

## # Hessian of log-density with respect to parameters.
## NormalToeplitz$set("public", "hess", overwrite = TRUE, function(z, dz, d2z, acf, dacf, d2acf) {
##   check_var(z, dx = dz, d2x = d2z, N = private$size, varname = "z")
##   check_var(acf, dx = dacf, d2x = d2acf, N = private$size, varname = "acf")
##   if(ncol(dz) != ncol(dacf)) {
##     stop("dz and dacf must have the same number of columns.")
##   }
##   if(!all(dim(d2z) == dim(d2acf))) {
##     stop("d2z and d2acf must have the same dimensions.")
##   }
##   ntheta <- ncol(dz)
##   .NormalToeplitz_hess(private$NTz_ptr, z, dz,
##                        matrix(d2z, private$size, ntheta*ntheta),
##                        acf, dacf,
##                        matrix(d2acf, private$size, ntheta*ntheta), ntheta)
## })

## # Full gradient of log-density.
## NormalToeplitz$set("public", "grad_full", overwrite = TRUE, function(z, acf) {
##   check_var(z, N = private$size, varname = "z")
##   check_var(acf, N = private$size, varname = "acf")
##   .NormalToeplitz_grad_full(private$NTz_ptr, z, acf)
## })

#--- formatting functions ------------------------------------------------------

# Check inputs to NormalToeplitz.
#
# @param x Main input variable.
# @param dx Optional gradient variable.
# @param d2x Optional Hessian variable.
# @param varname Name of variable.
#
# @details
# - If only `x` provided, checks length.
# - If `dx` provided, checks that its a matrix with `nrow(dx) = length(x)`.
# - If `d2x` is provided, check that `dim(d2x) = c(length(x), ncol(dx), ncol(dx))`.
check_var <- function(x, dx, d2x, N, varname) {
  if(!(is.vector(x) && is.numeric(x))) {
    stop(paste0(varname, " must be a numeric vector."))
  }
  if(length(x) != N) {
    stop(paste0(varname, "has wrong length."))
  }
  if(!missing(dx)) {
    if(!(is.matrix(dx) && is.numeric(dx))) {
      stop(paste0("d", varname, " must be a numeric matrix."))
    }
    if(nrow(dx) != N) {
      stop(paste0("d", varname, " has wrong number of rows."))
    }
    p <- ncol(dx)
    if(!missing(d2x)) {
      if(!(is.array(d2x) && (length(dim(d2x)) == 3) && is.numeric(d2x))) {
        stop(paste0("d2", varname, " must be a numeric 3D array."))
      }
      if(!all(dim(d2x) == c(N, p, p))) {
        stop(paste0("d2", varname, " has incompatible dimensions with d",
                    varname, " and ", varname, "."))
      }
    }
    if(!missing(d2x) && missing(dx)) {
      stop(paste0("Cannot provide d2", varname, " without d", varname, "."))
    }
  }
}

#--- scratch -------------------------------------------------------------------

## # exported constructor
## #' @rdname NormalToeplitz
## #'
## #' @param n Size of the Toeplitz matrix.
## #' @return A `NormalToeplitz` object.
## #' @export
## NormalToeplitz <- function(n) {
##   Nt <- .NormalToeplitz$new(n)
##   Nt
## }

## # now make methods display arguments
## .DollarNames.NormalToeplitz <- function(x, pattern) {
##   grep(pattern, getRefClass(class(x))$methods(), value = TRUE)
## }
