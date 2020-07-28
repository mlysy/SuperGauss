#' Multivariate normal with Toeplitz variance matrix.
#'
#' @description Provides methods for the Normal-Toeplitz (NTz) distribution defined as
#' ```
#' z ~ NTz(acf)   <=>   z ~ Normal(0, toeplitz(acf)),
#' ```
#' i.e., for a multivariate normal with mean zero and variance `Tz = toeplitz(acf)`.
#'
#' @export
NormalToeplitz <- R6Class(
  classname = "NormalToeplitz",

  private = list(

    NTz_ = NULL,
    N_ = NA,

    # deep clone method.
    # required to create new Xptr for NTz_ at C++ level
    deep_clone = function(name, value) {
      switch(name,
             NTz_ = {
               NTz_new <- NormalToeplitz_ctor(private$N_)
               ## if(self$has_acf()) {
               ##   .Toeplitz_set_acf(Tz_new, self$get_acf())
               ## }
               NTz_new
             },
             value)
    }


  ),

  public = list(

    #' @description Class constructor.
    #'
    #' @param N Size of the NTz random vector.
    #' @return A `NormalToeplitz` object.
    initialize = function(N) {
      private$NTz_ <- NormalToeplitz_ctor(N)
      private$N_ <- N
    },

    #' @description Get the size of the NTz random vector.
    #'
    #' @return Size of the NTz random vector.
    size = function() {
      private$N_
    },

    #' @description Log-density function.
    #'
    #' @param z Density argument. A vector of length `N` or an `N x n_obs` matrix where each column is an `N`-dimensional observation.
    #' @param acf A vector of length `N` containing the autocorrelation (i.e., first row/column) of the Toeplitz variance matrix.
    #' @return A scalar or vector of length `n_obs` containing the log-density of the NTz evaluated at its arguments.
    logdens = function(z, acf) {
      z <- as.matrix(z)
      check_ntz(z[,1], N = private$N_, varname = "z")
      check_ntz(acf, N = private$N_, varname = "acf")
      NormalToeplitz_logdens(private$NTz_, z, acf)
    },

    #' @description Gradient of the log-density with respect to parameters.
    #'
    #' @param z Density argument.  A vector of length `N`.
    #' @param dz An `N x n_theta` matrix containing the gradient `dz/dtheta`.
    #' @param acf A vector of length `N` containing the autocorrelation of the Toeplitz variance matrix.
    #' @param dacf An `N x n_theta` matrix containing the gradient `dacf/dtheta`.
    #' @param full_out If `TRUE`, returns the log-density as well (see 'Value').
    #' @return A vector of length `n_theta` containing the gradient of the NTz log-density with respect to `theta`, or a list with elements `ldens` and `grad` consisting of the log-density and the gradient vector.
    grad = function(z, dz, acf, dacf, full_out = FALSE) {
      check_ntz(z, dx = dz, N = private$N_, varname = "z")
      check_ntz(acf, dx = dacf, N = private$N_, varname = "acf")
      if(ncol(dz) != ncol(dacf)) {
        stop("dz and dacf must have the same number of columns.")
      }
      ## n_theta <- ncol(dz)
      NormalToeplitz_grad(private$NTz_, z, dz, acf, dacf, full_out)
    },

    #' @description Hessian of log-density with respect to parameters.
    #'
    #' @param z Density argument.  A vector of length `N`.
    #' @param dz An `N x n_theta` matrix containing the gradient `dz/dtheta`.
    #' @param d2z An `N x n_theta x n_theta` array containing the Hessian `d^2z/dtheta^2`.
    #' @param acf A vector of length `N` containing the autocorrelation of the Toeplitz variance matrix.
    #' @param dacf An `N x n_theta` matrix containing the gradient `dacf/dtheta`.
    #' @param d2acf An `N x n_theta x n_theta` array containing the Hessian `dacf^2/dtheta^2`.
    #' @param full_out If `TRUE`, returns the log-density and its gradient as well (see 'Value').
    #' @return An `n_theta x n_theta` matrix containing the Hessian of the NTz log-density with respect to `theta`, or a list with elements `ldens`, `grad`, and `hess` consisting of the log-density, its gradient (a vector of size `n_theta`), and the Hessian matrix, respectively.
    hess = function(z, dz, d2z, acf, dacf, d2acf, full_out = FALSE) {
      check_ntz(z, dx = dz, d2x = d2z, N = private$N_, varname = "z")
      check_ntz(acf, dx = dacf, d2x = d2acf, N = private$N_,
                varname = "acf")
      if(ncol(dz) != ncol(dacf)) {
        stop("dz and dacf must have the same number of columns.")
      }
      if(!all(dim(d2z) == dim(d2acf))) {
        stop("d2z and d2acf must have the same dimensions.")
      }
      n_theta <- ncol(dz)
      NormalToeplitz_hess(private$NTz_, z, dz,
                          matrix(d2z, private$N_, n_theta*n_theta),
                          acf, dacf,
                          matrix(d2acf, private$N_, n_theta*n_theta),
                          full_out)
    },

    #' @description Full gradient of log-density function.
    #'
    #' @param z Density argument.  A vector of length `N`.
    #' @param acf A vector of length `N` containing the autocorrelation of the Toeplitz variance matrix.
    #' @param calc_dldz Whether or not to calculate the gradient with respect to `z`.
    #' @param calc_dlda Whether or not to calculate the gradient with respect to `acf`.
    #' @return A list with elements:
    #' \describe{
    #'   \item{`ldens`}{The log-density evaluated at `z` and `acf`.}
    #'   \item{`dldz`}{The length-`N` gradient vector with respect to `z`, if `calc_dldz = TRUE`.}
    #'   \item{`dlda`}{The length-`N` gradient vector with respect to `acf`, if `calc_dlda = TRUE`.}
    #' }
    grad_full = function(z, acf, calc_dldz = TRUE, calc_dlda = TRUE) {
      if(!calc_dldz && !calc_dlda) {
        stop("At least one of calc_dldz or calc_dlda must be TRUE.")
      }
      check_ntz(z, N = private$N_, varname = "z")
      check_ntz(acf, N = private$N_, varname = "acf")
      NormalToeplitz_grad_full(private$NTz_, z, acf, calc_dldz, calc_dlda)
    }
  )

)

#--- formatting functions ------------------------------------------------------

#' Check inputs to NormalToeplitz.
#'
#' @param x Main input variable.
#' @param dx Optional gradient variable.
#' @param d2x Optional Hessian variable.
#' @param varname Name of variable.
#'
#' @details
#' - If only `x` provided, checks length.
#' - If `dx` provided, checks that its a matrix with `nrow(dx) = length(x)`.
#' - If `d2x` is provided, check that `dim(d2x) = c(length(x), ncol(dx), ncol(dx))`.
#' @noRd
check_ntz <- function(x, dx, d2x, N, varname) {
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

## # make methods display arguments
## # is this still needed?
## .DollarNames.NormalToeplitz <- function(x, pattern) {
##   grep(pattern, getRefClass(class(x))$methods(), value = TRUE)
## }

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

## # Class constructor
## NormalToeplitz$set("public", "initialize", overwrite = TRUE, function(n) {
##   private$NTz_ptr <- .NormalToeplitz_constructor(n)
##   private$size <- n
## })

## # Log-density function.
## NormalToeplitz$set("public", "logdens", overwrite = TRUE, function(z, acf) {
##   check_ntz(z, N = private$size, varname = "z")
##   check_ntz(acf, N = private$size, varname = "acf")
##   .NormalToeplitz_logdens(private$NTz_ptr, z, acf)
## })

## # Gradient of log-density with respect to parameters.
## NormalToeplitz$set("public", "grad", overwrite = TRUE, function(z, dz, acf, dacf) {
##   check_ntz(z, dx = dz, N = private$size, varname = "z")
##   check_ntz(acf, dx = dacf, N = private$size, varname = "acf")
##   if(ncol(dz) != ncol(dacf)) {
##     stop("dz and dacf must have the same number of columns.")
##   }
##   ntheta <- ncol(dz)
##   .NormalToeplitz_grad(private$NTz_ptr, z, dz, acf, dacf, ntheta)
## })

## # Hessian of log-density with respect to parameters.
## NormalToeplitz$set("public", "hess", overwrite = TRUE, function(z, dz, d2z, acf, dacf, d2acf) {
##   check_ntz(z, dx = dz, d2x = d2z, N = private$size, varname = "z")
##   check_ntz(acf, dx = dacf, d2x = d2acf, N = private$size, varname = "acf")
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
##   check_ntz(z, N = private$size, varname = "z")
##   check_ntz(acf, N = private$size, varname = "acf")
##   .NormalToeplitz_grad_full(private$NTz_ptr, z, acf)
## })
