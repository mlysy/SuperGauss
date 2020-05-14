#' Constructor and methods for Toeplitz matrix objects.
#'
#' @name Toeplitz
#'
#' @param x An R object.
#' @description The `Toeplitz` class contains efficient methods for linear algebra with symmetric positive definite (i.e., variance) Toeplitz matrices.
#'
#' @details An `N x N` Toeplitz matrix `Tz` is defined by its length-`N` "autocorrelation" vector `acf`, i.e., first row/column `Tz`.  Thus, for the function [stats::toeplitz()], we have `Tz = toeplitz(acf)`.
#'
#' It is assumed that `acf` defines a valid (i.e., positive definite) variance matrix.  The matrix multiplication methods still work when this is not the case but the other methods do not (return values typically contain `NaN`s).
#'
#' `as.Toeplitz(x)` attempts to convert its argument to a `Toeplitz` object by calling `Toeplitz$new(acf = x)`. `is.Toeplitz(x)` checks whether its argument is a `Toeplitz` object.
#'
#' @example examples/Toeplitz.R
NULL

#' @rdname Toeplitz
#' @export
is.Toeplitz <- function(x) "Toeplitz" %in% class(x)

#' @rdname Toeplitz
#' @export
as.Toeplitz <- function(x) {
  Toeplitz$new(acf = x)
}

#' @rdname Toeplitz
#' @export
Toeplitz <- R6Class(
  classname = "Toeplitz",

  private = list(

    Tz_ = NULL,
    PCG_ = NULL,
    N_ = NA,

    # deep clone method.
    # required to create new Xptr for Tz_ at C++ level
    deep_clone = function(name, value) {
      switch(name,
             Tz_ = {
               Tz_new <- .Toeplitz_constructor(private$N_)
               if(self$has_acf()) {
                 .Toeplitz_set_acf(Tz_new, self$get_acf())
               }
               Tz_new
             },
             PCG_ = .PCG_constructor(private$N_),
             value)
    }

  ),

  public = list(

    #' @description Class constructor.
    #'
    #' @param N Size of Toeplitz matrix.
    #' @param acf Autocorrelation vector of length `N`.
    #'
    #' @return A `Toeplitz` object.
    initialize = function(N, acf) {
      if(missing(N)) N <- length(acf)
      private$N_ <- N
      private$Tz_ <- .Toeplitz_constructor(N)
      private$PCG_ <- .PCG_constructor(N)
      if(!missing(acf)) self$set_acf(acf)
    },

    #' @description Print method.
    print = function() {
      if(self$has_acf()) {
        obj_acf <- self$get_acf()[1:min(6, private$N_)]
        obj_acf <- signif(obj_acf, digits = 3)
        if(private$N_ > 6) obj_acf <- c(obj_acf, "...")
      } else {
        obj_acf <- "NULL"
      }
      cat("Toeplitz matrix of size", private$N_, "\n",
          "acf: ", obj_acf, "\n")
    },

    #' @description Get the size of the Toeplitz matrix.
    #'
    #' @return Size of the Toeplitz matrix.  [`ncol()`][base::ncol()], [`nrow()`][base::nrow()], and [`dim()`][base::dim()] methods for `Toeplitz` objects also work as expected.
    size = function() {
      private$N_
    },

    #' @description Set the autocorrelation of the Toeplitz matrix.
    #'
    #' @param acf Autocorrelation vector of length `N`.
    set_acf = function(acf) {
      check_tz(acfs = list(acf = acf), N = private$N_)
      .Toeplitz_set_acf(private$Tz_, acf)
    },

    #' @description Get the autocorrelation of the Toeplitz matrix.
    #'
    #' @return The autocorrelation vector of length `N`.
    get_acf = function() {
      check_tz(has_acf = self$has_acf())
      .Toeplitz_get_acf(private$Tz_)
    },

    #' @description Check whether the autocorrelation of the Toeplitz matrix has been set.
    #'
    #' @return Logical; `TRUE` if `Toeplitz$set_acf()` has been called.
    has_acf = function() {
      .Toeplitz_has_acf(private$Tz_)
    },

    #' @description Toeplitz matrix-matrix product.
    #'
    #' @param x Vector or matrix with `N` rows.
    #' @return The matrix product `Tz %*% x`. `Tz %*% x` and `x %*% Tz` also work as expected.
    prod = function(x) {
      check_tz(has_acf = self$has_acf())
      if(is.vector(x)) x <- as.matrix(x)
      if(!(is.matrix(x) && is.numeric(x))) {
        stop("x must be a numeric matrix.")
      }
      if(nrow(x) != private$N_) {
        stop("Incompatible matrix multiplication dimensions.")
      }
      .Toeplitz_prod(private$Tz_, x)
    },

    #' @description Solve a Toeplitz system of equations.
    #'
    #' @param x Optional vector or matrix with `N` rows.
    #' @param method Solve method to use.  Choices are: `gschur` for a modified version of the Generalized Schur algorithm of Ammar & Gragg (1988), or `pcg` for the preconditioned conjugate gradient method of Chen et al (2006).  The former is faster and obtains the log-determinant as a direct biproduct.  The latter is more numerically stable for long-memory autocorrelations.
    #' @param tol Tolerance level for the `pcg` method.
    #' @return The solution in `z` to the system of equations `Tz %*% z = x`.  If `x` is missing, returns the inverse of `Tz`.  `solve(Tz, x)` and `solve(Tz, x, method, tol)` also work as expected.
    solve = function(x, method = c("gschur", "pcg"), tol = 1e-10) {
      method <- match.arg(method)
      check_tz(has_acf = self$has_acf())
      if(missing(x)) x <- diag(private$N_)
      vec_x <- is.vector(x)
      if(vec_x) x <- as.matrix(x)
      if(!(is.matrix(x) && is.numeric(x))) {
        stop("x must be a numeric matrix.")
      }
      if(nrow(x) != private$N_) {
        stop("Incompatible matrix solve dimensions.")
      }
      y <- switch(method,
                  gschur = .Toeplitz_solve(private$Tz_, x),
                  pcg = .PCG_solve(private$PCG_, self$get_acf(), x, tol))
      if(vec_x) y <- drop(y)
      y
    },

    #' @description Calculate the log-determinant of the Toeplitz matrix.
    #'
    #' @return The log-determinant `log(det(Tz))`.  `determinant(Tz)` also works as expected.
    log_det = function() {
      check_tz(has_acf = self$has_acf())
      .Toeplitz_log_det(private$Tz_)
    },

    #' @description Computes the trace-gradient with respect to Toeplitz matrices.
    #' @param acf2 Length-`N` autocorrelation vector of the second Toeplitz matrix.  This matrix must be symmetric but not necessarily positive definite.
    #' @return Computes the trace of
    #' ```
    #' solve(Tz, toeplitz(acf2)).
    #' ```
    #' This is used in the computation of the gradient of `log(det(Tz(theta)))` with respect to `theta`.
    trace_grad = function(acf2) {
      check_tz(acfs = list(acf2 = acf2), N = private$N_,
               has_acf = self$has_acf())
      .Toeplitz_trace_grad(private$Tz_, acf2)
    },

    #' @description Computes the trace-Hessian with respect to Toeplitz matrices.
    #'
    #' @param acf2 Length-`N` autocorrelation vector of the second Toeplitz matrix.  This matrix must be symmetric but not necessarily positive definite.
    #' @param acf3 Length-`N` autocorrelation vector of the third Toeplitz matrix.  This matrix must be symmetric but not necessarily positive definite.
    #' @return Computes the trace of
    #' ```
    #' solve(Tz, toeplitz(acf2)) %*% solve(Tz, toeplitz(acf3)).
    #' ```
    #' This is used in the computation of the Hessian of `log(det(Tz(theta)))` with respect to `theta`.
    trace_hess = function(acf2, acf3) {
      check_tz(acfs = list(acf2 = acf2, acf3 = acf3),
               N = private$N_, has_acf = self$has_acf())
      .Toeplitz_trace_hess(private$Tz_, acf2, acf3)
    }

  )
)
setOldClass("Toeplitz")

#--- generic methods -----------------------------------------------------------

# ncol
#' @rdname Toeplitz
#' @aliases ncol,Toeplitz-method
#' @usage NULL
#' @export
setMethod("ncol", "Toeplitz", function(x) {
  x$size()
})

# nrow
#' @rdname Toeplitz
#' @aliases nrow,Toeplitz-method
#' @usage NULL
#' @export
setMethod("nrow", "Toeplitz", function(x) {
  x$size()
})

# dim
#' @rdname Toeplitz
#' @export
dim.Toeplitz <- function(x) rep(x$size(), 2)

# Matrix multiplication
#' @rdname Toeplitz
#' @aliases %*%
#' @usage NULL
#' @export
`%*%` <- function(x, y) UseMethod("%*%")

#' @export
`%*%.default` <- function(x, y) {
  if(is.Toeplitz(x) || is.Toeplitz(y)) return(`%*%.Toeplitz`(x, y))
  base::`%*%`(x, y)
}

#' @export
`%*%.Toeplitz` <- function(x, y) {
  if(is.Toeplitz(x)) {
    # Toeplitz %*% Matrix
    if(is.Toeplitz(y)) {
      # Toeplitz %*% Toeplitz
      # This can be done more efficiently by Gohberg-Semencul decomposition,
      # but not currently implemented.
      y <- stats::toeplitz(y$get_acf())
    }
    check_tz(has_acf = x$has_acf())
    if(is.vector(y)) y <- as.matrix(y)
    if(!(is.matrix(y) && is.numeric(y))) {
      stop("Second argument must be a numeric matrix.")
    }
    if(nrow(y) != x$size()) {
      stop("Incompatible matrix multiplication dimensions.")
    }
    ans <- x$prod(y)
  } else if(is.Toeplitz(y)) {
    # Matrix %*% Toeplitz
    check_tz(has_acf = y$has_acf())
    if(is.vector(x)) x <- as.matrix(x)
    if(!(is.matrix(x) && is.numeric(x))) {
      stop("First argument must be a numeric matrix.")
    }
    if(ncol(x) != y$size()) {
      stop("Incompatible matrix multiplication dimensions.")
    }
    ans <- t(y$prod(t(x)))
  } else stop("This shouldn't happen.  Please contact package maintainer.")
  ans
}

# determinant
#' @rdname Toeplitz
#' @aliases determinant determinant,Toeplitz-method
#' @usage NULL
#' @export
setMethod("determinant", "Toeplitz", function(x, logarithm = TRUE, ...) {
  ldT <- x$log_det()
  if(!logarithm) ldT <- exp(ldT)
  ldT
})

# solve
#' @rdname Toeplitz
#' @aliases solve solve,Toeplitz,ANY-method solve,Toeplitz-method
#' @usage NULL
#' @export
setMethod(
  f = "solve",
  signature = "Toeplitz",
  definition = function(a, b, method = c("gschur", "pcg"), tol = 1e-10, ...) {
    check_tz(has_acf = a$has_acf())
    if(missing(b)) b <- diag(a$size())
    vec_b <- is.vector(b)
    if(vec_b) b <- as.matrix(b)
    if(!is.matrix(b)) {
      stop("b must be a matrix or vector.")
    }
    if(nrow(b) != a$size()) {
      stop("Incompatible matrix solve dimensions.")
    }
    y <- a$solve(b, method, tol)
    if(vec_b) y <- drop(y)
    y
  }
)


#--- helper functions ----------------------------------------------------------

#' Check inputs to Toeplitz.
#'
#' @param acfs Named list of acfs.
#' @param N Required size of each acf.
#' @param has_acf Optional whether or not acf has been set.
#'
#' @details
#' - For each element of `acfs` checks that it has length `N`.
#' - If `has_acf` provided returns error if `has_acf == FALSE`.
#' @noRd
check_tz <- function(acfs, N, has_acf) {
  if(!missing(has_acf)) {
    # check has_acf
    if(!has_acf) stop("Toeplitz$set_acf() has not been called yet.")
  }
  if(!missing(acfs)) {
    # check acfs
    for(varname in names(acfs)) {
      x <- acfs[[varname]]
      if(!(is.vector(x) && is.numeric(x))) {
        stop(paste0(varname, " must be a numeric vector."))
      }
      if(length(x) != N) {
        stop(paste0(varname, "has wrong length."))
      }
    }
  }
}

