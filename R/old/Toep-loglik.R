#' Log-likelihood function for LMN models using Superfast Algorithm
#'
#' @export
toep.loglik <- function(Beta, Sigma, suff, Y, X, V, acf, Toep = Toep) {
  if(missing(suff)) {
    suff <- toep.suff(Y = Y, X = X, V = V, acf = acf, Toep = Toep)
  }
  # get sufficient statistics for likelihood evaluation
  n <- suff$n
  S <- suff$S
  q <- nrow(S)
  ldV <- suff$ldV
  Beta.hat <- suff$Beta.hat
  T <- suff$T
  noBeta <- is.null(Beta.hat)
  # log-likelihood calculation
  if(!noBeta) {
    Z <- Beta-Beta.hat
    S <- S + crossprod(Z, T %*% Z)
  }
  -0.5 * (sum(diag(solve(Sigma,S))) + q*ldV + n*ldet(Sigma) + n*q * log(2*pi))
}
