#' Constructor and methods for Toeplitz matrix objects.
#'
#' The `Toeplitz` class contains efficient methods for linear algebra with symmetric positive definite (i.e., variance) Toeplitz matrices.
#'
#' @aliases Toeplitz set_acf get_acf trace_grad trace_hess show.Toeplitz %*% determinant solve %*%,ANY,Toeplitz-method %*%,Toeplitz,ANY-method determinant,Toeplitz-method dim.Toeplitz ncol,Toeplitz-method nrow,Toeplitz-method show,Toeplitz-method solve,Toeplitz-method
#' @section Methods:
#' If `Tz` is a `Toeplitz` object with first row/column given by `acf`, then:
#' \describe{
#' \item{`Tz$set_acf(acf)`}{Sets the autocorrelation of the matrix.}
#' \item{`Tz$get_acf()`}{Gets the autocorrelation of the matrix.}
#' \item{`nrow(Tz)`, `ncol(Tz)`, `dim(Tz)`}{Selected dimension(s) of the matrix.}
#' \item{`Tz %*% X`, `X %*% Tz`}{Toeplitz-Matrix and Matrix-Toeplitz multiplication.  Also works if `X` is a vector.}
#' \item{`solve(Tz, X)`, `solve(Tz)`}{Solves Toeplitz systems of equations.  When second argument is missing, returns the inverse of the Toeplitz matrix.}
#' \item{`determinant(Tz)`}{Log-determinant of the Toeplitz matrix, i.e., same thing as `log(det(toeplitz(acf)))`.}
#' \item{`Tz$trace_grad(acf2)`}{Computes the trace of `solve(Tz, toeplitz(acf2))`.  This is used in the computation of the gradient of `log(det(Tz(theta)))` with respect to `theta`.}
#' \item{`Toep$trace_hess(acf2, acf3)`}{Computes the trace of
#' ```
#' solve(Tz, toeplitz(acf2)) %*% solve(Tz, toeplitz(acf3)).
#' ```
#' This is used in the computation of the Hessian of `log(det(Tz(theta)))` with respect to `theta`.}
#' }
#'
#' \subsection{`Tz$set_acf(acf)`}{
#'
#' Sets the autocorrelation of the matrix.
#'
#' }
#'
#' \subsection{`Tz$trace_hess(acf2, acf3)`}{
#'
#' Computes the trace of
#' ```
#' solve(Tz, toeplitz(acf2)) %*% solve(Tz, toeplitz(acf3)).
#' ```
#' This is used in the computation of the Hessian of `log(det(Tz(theta)))` with respect to `theta`.
#' }
#'
#' @details It is assumed that the autocorrelation of the `Toeplitz` object defines a valid (i.e., positive definite) variance matrix.  The multiplication algorithms still work when this is not the case but the other algorithms do not (return values typically contain `NaN`s).
#' @examples
#' # construction
#' acf <- exp(-(1:5))
#' Toep <- Toeplitz(acf = acf)
#' # alternatively, can allocate space first
#' Toep <- Toeplitz(n = length(acf))
#' Toep$set_acf(acf = acf)
#'
#' dim(Toep) # == c(nrow(Toep), ncol(Toep))
#' Toep # show method
#' Toep$get_acf() # extract the acf
#'
#' # linear algebra
#' X <- matrix(rnorm(10), 5, 2)
#' Toep %*% X
#' t(X) %*% Toep
#' solve(Toep, X)
#' determinant(Toep) # log-determinant
#' @export
Toeplitz <- R6Class(
  classname = "Toeplitz",

  private = list(

    Tz_ptr = NULL,
    N_ = NA

  ),

  active = list(
    has_acf = function() {
      .Toeplitz_has_acf(private$Tz_ptr)
    }
  ),

  public = list(

    initialize = function(N) {
      private$N_ <- N
      private$Tz_ptr <- .Toeplitz_constructor(N)
    },

    print = function(...) {
      if(self$has_acf) {
        obj_acf <- self$get_acf()[1:min(6, private$N_)]
        obj_acf <- signif(obj_acf, digits = 3)
        if(private$N_ > 6) obj_acf <- c(obj_acf, "...")
      } else {
        obj_acf <- "NULL"
      }
      cat("Toeplitz matrix of size", private$N_, "\n",
          "acf: ", obj.acf, "\n")
    },

    size = function() {
      private$N_
    },

    set_acf = function(acf) {
      check_tz(acfs = list(acf = acf), N = private$N_)
      .Toeplitz_set_acf(private$Tz_ptr, acf)
    },

    get_acf = function() {
      check_tz(has_acf = self$has_acf)
      .Toeplitz_get_acf(private$Tz_ptr)
    },

    prod = function(x) {
      check_tz(has_acf = self$has_acf)
      if(is.vector(x)) x <- as.matrix(x)
      if(!(is.matrix(x) && is.numeric(x))) {
        stop("x must be a numeric matrix.")
      }
      if(nrow(x) != private$N_) {
        stop("Incompatible matrix multiplication dimensions.")
      }
      .Toeplitz_prod(private$Tz_ptr, x)
    },

    solve = function(x) {
      check_tz(has_acf = self$has_acf)
      if(missing(x)) x <- diag(private$N_)
      if(is.vector(x)) x <- as.matrix(x)
      if(!(is.matrix(x) && is.numeric(x))) {
        stop("x must be a numeric matrix.")
      }
      if(nrow(x) != private$N_) {
        stop("Incompatible matrix solve dimensions.")
      }
      .Toeplitz_solve(private$Tz_ptr, x)
    },

    log_det = function() {
      check_tz(has_acf = self$has_acf)
      .Toeplitz_log_det(private$Tz_ptr)
    },

    trace_grad = function(acf2) {
      check_tz(acfs = list(acf2 = acf2), N = private$N_,
               has_acf = self$has_acf)
      .Toeplitz_trace_grad(private$Tz_ptr, acf2)
    },

    trace_hess = function(acf2, acf3) {
      check_tz(acfs = list(acf2 = acf2, acf3 = acf3),
               N = private$N_, has_acf = self$has_acf)
      .Toeplitz_trace_hess(private$Tz_ptr, acf2, acf3)
    }

  )
)
setOldClass("Toeplitz")

#--- generic methods -----------------------------------------------------------

## # show
## #' @export
## setMethod("show", "Toeplitz", function(object) {
##   if(object$has_acf) {
##     obj.acf <- object$get_acf()[1:min(6, object$size)]
##     obj.acf <- signif(obj.acf, digits = 3)
##     if(object$size() > 6) obj.acf <- c(obj.acf, "...")
##   } else {
##     obj.acf <- "NULL"
##   }
##   cat("Toeplitz matrix of size", object$size(), "\n",
##       "acf: ", obj.acf, "\n")
## })

is.Toeplitz <- function(x) "Toeplitz" %in% class(x)

# ncol
#' @export
setMethod("ncol", "Toeplitz", function(x) {
  x$size()
})

# nrow
#' @export
setMethod("nrow", "Toeplitz", function(x) {
  x$size()
})

# dim
#' @export
dim.Toeplitz <- function(x) rep(x$size(), 2)
## setMethod("dim", "Toeplitz", function(x) {
##   c(x$size(), x$size())
## })

# Matrix multiplication
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
    if(is.Toeplitz(y)) y <- toeplitz(y$get_acf())
    check_tz(has_acf = x$has_acf)
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
    check_tz(has_acf = y$has_acf)
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

## # Toeplitz-Matrix multiplication
## #' @export
## setMethod("%*%", signature(x = "Toeplitz", y = "ANY"), function(x, y) {
##   check_tz(has_acf = x$has_acf)
##   if(is.vector(y)) y <- as.matrix(y)
##   if(!(is.matrix(y) && is.numeric(y))) {
##     stop("Second argument must be a numeric matrix.")
##   }
##   if(nrow(y) != x$size()) {
##     stop("Incompatible matrix multiplication dimensions.")
##   }
##   x$prod(y)
## })

## # Matrix-Toeplitz multiplication
## #' @export
## setMethod("%*%", signature(x = "ANY", y = "Toeplitz"), function(x, y) {
##   check_tz(has_acf = x$has_acf)
##   if(is.vector(x)) x <- as.matrix(x)
##   if(!(is.matrix(x) && is.numeric(x))) {
##     stop("First argument must be a numeric matrix.")
##   }
##   if(ncol(x) != y$size()) {
##     stop("Incompatible matrix multiplication dimensions.")
##   }
##   t(y$prod(t(x)))
## })

# determinant
#' @export
setMethod("determinant", "Toeplitz", function(x, logarithm = TRUE, ...) {
  ## if(!.Toeplitz_has_acf(x$cpp_ptr)) {
  ##   stop("set_acf has not been called yet")
  ## }
  ## ldT <- .Toeplitz_log_det(x$cpp_ptr)
  ldT <- x$log_det()
  if(!logarithm) ldT <- exp(ldT)
  ldT
})

# solve
#' @export
setMethod("solve", "Toeplitz", function(a, b, ...) {
  check_tz(has_acf = a$has_acf)
  ## if(!.Toeplitz_has_acf(a$cpp_ptr)) {
  ##   stop("set_acf has not been called yet")
  ## }
  if(missing(b)) b <- diag(a$size())
  if(is.vector(b)) b <- as.matrix(b)
  if(!is.matrix(b)) {
    stop("b must be a matrix or vector.")
  }
  if(nrow(b) != a$size()) {
    stop("Incompatible matrix solve dimensions.")
    ## stop("a and b are non-conformable.")
  }
  a$solve(b)
  ## .Toeplitz_solve(a$cpp_ptr, b)
})


#--- helper functions ----------------------------------------------------------

#' Check inputs to Toeplitz.
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

#--- old -----------------------------------------------------------------------

## if(FALSE) {
## .Toeplitz <- setRefClass("Toeplitz",
##                          fields = list(cpp_ptr = "externalptr",
##                                        size = "numeric"))
## .Toeplitz$lock("cpp_ptr") # locked fields
## .Toeplitz$lock("size")
## # internal constructor
## .Toeplitz$methods(initialize = function(n) {
##   cpp_ptr <<- .Toeplitz_constructor(n)
##   size <<- n
## })

## # exported constructor
## #' @rdname Toeplitz-class
## #' @usage Toeplitz(n, acf)
## #' @usage Toeplitz$set_acf(acf)
## #' @usage Toeplitz$get_acf()
## #' @usage Toeplitz$trace_grad(acf2)
## #' @usage Toeplitz$trace_hess(acf2, acf3)
## #' @param n Size of the Toeplitz matrix.
## #' @param acf Autocorrelation vector of Toeplitz matrix.
## #' @param acf2 Autocorrelation of second Toeplitz matrix.
## #' @param acf3 Autocorrelation of third Toeplitz matrix.
## #' @return A `Toeplitz` object.
## #' @export
## Toeplitz <- function(n, acf) {
##   if(missing(n)){
##     n <- length(acf)
##   }
##   Tz <- .Toeplitz$new(n)
##   if(!missing(acf)) {
##     Tz$set_acf(acf)
##   }
##   Tz
## }

## #--- custom methods ------------------------------------------------------------

## # setter
## .Toeplitz$methods(set_acf = function(acf) {
##   if(length(acf) != size) {
##     stop("acf has wrong length.")
##   }
##   .Toeplitz_set_acf(cpp_ptr, acf)
## })

## # getter
## .Toeplitz$methods(get_acf = function() {
##   .Toeplitz_get_acf(cpp_ptr)
## })

## # trace_grad
## .Toeplitz$methods(trace_grad = function(acf2) {
##   if(!.Toeplitz_has_acf(cpp_ptr)) {
##     stop("set_acf has not been called yet")
##   }
##   if(length(acf2) != size) {
##     stop("acf2 has wrong length.")
##   }
##   .Toeplitz_trace_grad(cpp_ptr, acf2)
## })

## # trace_hess
## .Toeplitz$methods(trace_hess = function(acf2, acf3) {
##   if(!.Toeplitz_has_acf(cpp_ptr)) {
##     stop("set_acf has not been called yet")
##   }
##   if(length(acf2) != size) {
##     stop("acf2 has wrong length.")
##   }
##   if(length(acf3) != size) {
##     stop("acf3 has wrong length.")
##   }
##   .Toeplitz_trace_hess(cpp_ptr, acf2, acf3)
## })


## #--- generic methods -----------------------------------------------------------

## # show
## #' @export
## setMethod("show", "Toeplitz", function(object) {
##   if(.Toeplitz_has_acf(object$cpp_ptr)) {
##     obj.acf <- object$get_acf()[1:min(6, object$size)]
##     obj.acf <- signif(obj.acf, digits = 3)
##     if(object$size > 6) obj.acf <- c(obj.acf, "...")
##   } else {
##     obj.acf <- "NULL"
##   }
##   cat("Toeplitz matrix of size", object$size, "\n",
##       "acf: ", obj.acf, "\n")
## })

## # ncol
## #' @export
## setMethod("ncol", "Toeplitz", function(x) {
##   x$size
## })

## # nrow
## #' @export
## setMethod("nrow", "Toeplitz", function(x) {
##   x$size
## })

## # dim
## #' @export
## setMethod("dim", "Toeplitz", function(x) {
##   rep(x$size, 2)
## })

## # Toeplitz-Matrix multiplication
## #' @export
## setMethod("%*%", signature(x = "Toeplitz", y = "ANY"), function(x, y) {
##   if(!.Toeplitz_has_acf(x$cpp_ptr)) {
##     stop("set_acf has not been called yet")
##   }
##   if(is.vector(y)) y <- as.matrix(y)
##   if(!is.matrix(y)) {
##     stop("Second argument should be a matrix or vector.")
##   }
##   if(nrow(y) != x$size) {
##     stop("Toeplitz and second argument are non-conformable.")
##   }
##   .Toeplitz_prod(x$cpp_ptr, y)
## })

## # Matrix-Toeplitz multiplication
## setMethod("%*%", signature(x = "ANY", y = "Toeplitz"), function(x, y) {
##   if(!.Toeplitz_has_acf(y$cpp_ptr)) {
##     stop("set_acf has not been called yet")
##   }
##   if(is.vector(x)) x <- as.matrix(x)
##   if(!is.matrix(x)) {
##     stop("First argument should be a matrix or vector.")
##   }
##   if(ncol(x) != y$size) {
##     stop("First argument and Toeplitz are non-conformable.")
##   }
##   t(.Toeplitz_prod(y$cpp_ptr, t(x)))
## })

## # determinant
## #' @export
## setMethod("determinant", "Toeplitz",
##           function(x, logarithm = TRUE, ...) {
##   if(!.Toeplitz_has_acf(x$cpp_ptr)) {
##     stop("set_acf has not been called yet")
##   }
##   ldT <- .Toeplitz_log_det(x$cpp_ptr)
##   if(!logarithm) {
##     ldT <- exp(ldT)
##   }
##   ldT
## })

## # solve
## #' @export
## setMethod("solve", "Toeplitz", function(a, b, ...) {
##   if(!.Toeplitz_has_acf(a$cpp_ptr)) {
##     stop("set_acf has not been called yet")
##   }
##   if(missing(b)) b <- diag(a$size)
##   if(is.vector(b)) b <- as.matrix(b)
##   if(!is.matrix(b)) {
##     stop("b must be a matrix or vector.")
##   }
##   if(nrow(b) != a$size) {
##     stop("a and b are non-conformable.")
##   }
##   .Toeplitz_solve(a$cpp_ptr, b)
## })

## # now make methods display arguments
## .DollarNames.Toeplitz <- function(x, pattern)
##     grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
## }

## #--- new method ----------------------------------------------------------------


## ## #' @exportClass Toeplitz_cpp
## ## Rcpp::setRcppClass("Toeplitz_cpp", "Toeplitz_cpp", module = "Toeplitz")

## ## #' new Toeplitz Class
## ## #'
## ## #' @export
## ## .Toeplitz_R <- setRefClass("Toeplitz_R",
## ##                            contains = "Toeplitz_cpp",
## ##                            fields = list(size = "integer"),
## ##                            methods = list(
## ##                              initialize = function(n) {
## ##                                callSuper(n)
## ##                                size <<- n
## ##                              },
## ##                              set_acf = function(acf) {
## ##                                if(length(acf) != size) {
## ##                                  stop("acf has wrong length.")
## ##                                }
## ##                                callSuper(acf)
## ##                              }
## ##                            ))

## ## #' New exported constructor
## ## #'
## ## #' @export
## ## Toeplitz_R <- function(n, acf) {
## ##   if(missing(n)){
## ##     n <- length(acf)
## ##   }
## ##   Tz <- .Toeplitz_R$new(n)
## ##   if(!missing(acf)) {
## ##     Tz$set_acf(acf)
## ##   }
## ##   Tz
## ## }


