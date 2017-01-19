#' Wrapper for the log-determinant of a matrix
#' @param V A square matrix.
#' @return The log-determinant \code{log(det(V))}.
#' @keywords internal
ldet <- function(V) {
  determinant(V, log = TRUE)$mod[1]
}

# Log of Multi-Gamma Function
lmgamma <- function(x, p) {
  p*(p-1)/4 * log(pi) + sum(lgamma(x + (1-1:p)/2))
}

# Non-Symmetric Toeplitz Matrix
toeplitz2 <- function(col, row, debug = FALSE) {
  # dimensions
  n <- length(col)
  d <- length(row)
  if(col[1] != row[1]) {
    stop("row[1] and col[1] must be the same.")
  }
  CR <- c(rev(row[-1]), col)
  T <- matrix(NA, n, d)
  if(debug) browser()
  for(ii in 1:d) {
    T[,ii] <- CR[(d-ii)+1:n]
  }
  T
}

# names of prior
.PriorNames <- c("Lambda", "Omega", "Psi", "nu")

# Default prior specification
.DefaultPrior <- function(prior, p, q, noSigma) {
  noBeta <- p == 0
  # extract elements
  Lambda <- prior$Lambda
  Omega <- prior$Omega
  nu <- prior$nu
  Psi <- prior$Psi
  # assign defaults
  if(!noBeta) {
    noBeta <- (length(Omega) == 1) && is.na(Omega)
  }
  if(noBeta) {
    Lambda <- 0
    Omega <- NA
  } else {
    if(is.null(Lambda)) Lambda <- matrix(0,p,q)
    if(is.null(Omega) || all(Omega == 0)) Omega <- matrix(0,p,p)
  }
  if(missing(noSigma) || !noSigma) {
    noSigma <- (length(nu) == 1) && is.na(nu)
  }
  if(noSigma) {
    Psi <- 0
    nu <- NA
  } else {
    if(is.null(nu)) nu <- 0
    if(is.null(Psi) || all(Psi == 0)) Psi <- matrix(0,q,q)
  }
  list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
}
