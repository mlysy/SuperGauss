#' Multivariate normal with Circulant variance matrix.
#'
#' @description Provides methods for the Normal-Circulant (NCt) distribution, which for a random vector `z` of length `N` is defined as
#' ```
#' z ~ NCt(uacf)   <=>   z ~ Normal(0, toeplitz(acf)),
#' ```
#' where `uacf` are the `Nu = floor(N/2)+1` unique elements of the autocorrelation vector `acf`, i.e.,
#' ```
#' acf = (uacf, rev(uacf[2:(Nu-1)]),   N even,
#'     = (uacf, rev(uacf[2:Nu])),      N odd.
#' ```
#' @export
NormalCirculant <- R6Class(
  classname = "NormalCirculant",

  private = list(

    NCt_ = NULL,
    N_ = NA,
    Nu_ = NA,

    # deep clone method.
    # required to create new Xptr for NCt_ at C++ level
    deep_clone = function(name, value) {
      switch(name,
             NCt_ = {
               NCt_new <- NormalCirculant_ctor(private$N_)
               NCt_new
             },
             value)
    }


  ),

  public = list(

    #' @description Class constructor.
    #'
    #' @param N Size of the NCt random vector.
    #' @return A `NormalCirculant` object.
    initialize = function(N) {
      private$NCt_ <- NormalCirculant_ctor(N)
      private$N_ <- N
      private$Nu_ <- N%/%2 + 1
    },

    #' @description Get the size of the NCt random vector.
    #'
    #' @return Size of the NCt random vector.
    size = function() {
      private$N_
    },

    #' @description Log-density function.
    #'
    #' @param z Density argument. A vector of length `N` or an `N x n_obs` matrix where each column is an `N`-dimensional observation.
    #' @param uacf A vector of length `Nu = floor(N/2)` containing the first half of the autocorrelation (i.e., first row/column) of the Circulant variance matrix.
    #' @return A scalar or vector of length `n_obs` containing the log-density of the NCt evaluated at its arguments.
    logdens = function(z, uacf) {
      check_ntz(z, N = private$N_, varname = "z")
      check_ntz(uacf, N = private$Nu_, varname = "uacf")
      NormalCirculant_logdens(private$NCt_, z, uacf)
    },

    #' @description Full gradient of log-density function.
    #'
    #' @param z Density argument.  A vector of length `N`.
    #' @param uacf A vector of length `Nu = floor(N/2)` containing the first half of the autocorrelation (i.e., first row/column) of the Circulant variance matrix.
    #' @param calc_dldz Whether or not to calculate the gradient with respect to `z`.
    #' @param calc_dldu Whether or not to calculate the gradient with respect to `uacf`.
    #' @return A list with one or both elements `dldz` and `dldu`, corresponding to the length `N` gradient vector with respect to `z` and the length `Nu = floor(N/2)+1` gradient vector with respect to `uacf`.
    grad_full = function(z, uacf, calc_dldz = TRUE, calc_dldu = TRUE) {
      if(!calc_dldz && !calc_dldu) {
        stop("At least one of calc_dldz or calc_dldu must be TRUE.")
      }
      check_ntz(z, N = private$N_, varname = "z")
      check_ntz(uacf, N = private$Nu_, varname = "uacf")
      NormalCirculant_grad_full(private$NCt_, z, uacf,
                                calc_dldz, calc_dldu)
    }

  )

)

