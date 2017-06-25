# Wrapper for the log-determinant of a matrix
# @param V A square matrix.
# @return The log-determinant \code{log(det(V))}.
# @keywords internal
ldet <- function(V) {
  determinant(V, logarithm = TRUE)$mod[1]
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
