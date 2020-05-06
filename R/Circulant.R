#' Constructor and methods for Circulant matrix objects.
#'
#' @name Circulant
NULL

#' @rdname Circulant
#' @export
Circulant <- R6Class(
  classname = "Circulant",

  private = list(

    Ct_ = NULL,
    N_ = NA

  ),

  public = list(

    #' @description Class constructor.
    #'
    #' @param N Size of Circulant matrix.
    #' @param uacf Optional vector of `Nu = floor(N/2)+1` unique elements of the autocorrelation.
    #' @param upsd Optional vector of `Nu = floor(N/2)+1` unique elements of the PSD.
    #'
    #' @return A `Circulant` object.
    initialize = function(N, uacf, upsd) {
      private$N_ <- N
      private$Ct_ <- Circulant_ctor(private$N_)
      if(!missing(uacf)) {
        if(!missing(upsd)) stop("Cannot provide both uacf and upsd.")
        self$set_acf(uacf)
      }
      if(!missing(upsd)) {
        self$set_psd(upsd)
      }
    },

    #' @description Get the size of the Circulant matrix.
    #'
    #' @return Size of the Circulant matrix.
    size = function() {
      private$N_
    },

    #' @description Set the autocorrelation of the Circulant matrix.
    #'
    #' @param uacf Vector of `Nu = floor(N/2)+1` unique elements of the autocorrelation.
    set_acf = function(uacf) {
      check_ct(uvecs = list(uacf = uacf), N = private$N_)
      Circulant_set_acf(private$Ct_, uacf)
    },

    #' @description Get the autocorrelation of the Circulant matrix.
    #'
    #' @return The complete autocorrelation vector of length `N`.
    get_acf = function() {
      check_ct(has_acf = self$has_acf())
      Circulant_get_acf(private$Ct_)
    },

    #' @description Set the PSD of the Circulant matrix.
    #'
    #' The power spectral density (PSD) of a Circulant matrix `Ct = Circulant(acf)` is defined as `psd = iFFT(acf)`.
    #'
    #' @param upsd Vector of `Nu = floor(N/2)+1` unique elements of the psd.
    set_psd = function(upsd) {
      check_ct(uvecs = list(upsd = upsd), N = private$N_)
      Circulant_set_psd(private$Ct_, upsd)
    },

    #' @description Get the PSD of the Circulant matrix.
    #'
    #' @return The complete PSD vector of length `N`.
    get_psd = function() {
      check_ct(has_acf = self$has_acf())
      Circulant_get_psd(private$Ct_)
    },

    #' @description Check whether the autocorrelation of the Circulant matrix has been set.
    #'
    #' @return Logical; `TRUE` if `Circulant$set_acf()` has been called.
    has_acf = function() {
      Circulant_has_acf(private$Ct_)
    },

    #' @description Circulant matrix-matrix product.
    #'
    #' @param x Vector or matrix with `N` rows.
    #' @return The matrix product `Ct %*% x`.
    prod = function(x) {
      check_ct(has_acf = self$has_acf())
      if(is.vector(x)) x <- as.matrix(x)
      if(!(is.matrix(x) && is.numeric(x))) {
        stop("x must be a numeric matrix.")
      }
      if(nrow(x) != private$N_) {
        stop("Incompatible matrix multiplication dimensions.")
      }
      Circulant_prod(private$Ct_, x)
    },

    #' @description Solve a Circulant system of equations.
    #'
    #' @param x Optional vector or matrix with `N` rows.
    #' @return The solution in `z` to the system of equations `Ct %*% z = x`.  If `x` is missing, returns the inverse of `Ct`.
    solve = function(x) {
      check_ct(has_acf = self$has_acf())
      if(missing(x)) x <- diag(private$N_)
      vec_x <- is.vector(x)
      if(vec_x) x <- as.matrix(x)
      if(!(is.matrix(x) && is.numeric(x))) {
        stop("x must be a numeric matrix.")
      }
      if(nrow(x) != private$N_) {
        stop("Incompatible matrix solve dimensions.")
      }
      y <- Circulant_solve(private$Ct_, x)
      if(vec_x) y <- drop(y)
      y
    },

    #' @description Calculate the log-determinant of the Circulant matrix.
    #'
    #' @return The log-determinant `log(det(Ct))`.
    log_det = function() {
      check_ct(has_acf = self$has_acf())
      Circulant_log_det(private$Ct_)
    }

  )

)
setOldClass("Circulant")

#--- helper functions ----------------------------------------------------------

#' Check inputs to Circulant.
#'
#' @param uvecs Named list of unique-element vectors (usually uacf or upsd).
#' @param N Required size of Circulant matrix.
#' @param has_acf Optional whether or not acf has been set.
#'
#' @details
#' - For each element of `uvecs` checks that it has length `floor(N/2)+1`.
#' - If `has_acf` provided returns error if `has_acf == FALSE`.
#' @noRd
check_ct <- function(uvecs, N, has_acf) {
  if(!missing(has_acf)) {
    # check has_acf
    if(!has_acf) stop("Circulant$set_acf() has not been called yet.")
  }
  if(!missing(uvecs)) {
    # check uvecs
    for(varname in names(uvecs)) {
      x <- uvecs[[varname]]
      varname <- "uacf"
      if(!(is.vector(x) && is.numeric(x))) {
        stop(paste0(varname, " must be a numeric vector."))
      }
      Nu <- N%/%2 + 1
      if(length(x) != Nu) {
        stop(paste0(varname, "has wrong length."))
      }
    }
  }
}
