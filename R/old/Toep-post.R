#' Parameters of the Posterior Conditional MNIW Distribution
#'
#' @param prior A list with elements \code{Lambda, Omega, Psi, nu} of the prior MNIW distribution.  Any omitted or null value defaults to zero  (see details).
#' @param noSigma When \code{True} assumes that \code{Sigma = diag(q)}.
#' @param calc.prior When \code{True} also returns the parameters of the prior distribution.
#' @return A list with elements \code{Lambda, Omega, Psi, nu}, the parameters of the MNIW posterior distribution.  Also returns element \code{prior}, a list with the same elements for the conjugate prior distribution.
#'
#' \code{(Beta, Sigma) | Y ~ MNIW(Lambda, Omega, Psi, nu)}.
#'
#' Also returns \code{prior}, a list of the parameters of the MNIW prior.
#' @details The MNIW distribution is parametrized as
#' \deqn{Sigma ~ iWish(Psi, nu), \qquad Beta ~ MNorm(Lambda, Omega^{-1}, Sigma.}
#' As this is a posterior distribution, the parameter always define a proper MNIW distribution.  \code{nu = NA} and \code{Omega = NA} respectively mean that \code{Beta = 0} or \code{Sigma = diag(q)} is known.
#' @export
toep.post <- function(suff, Y, X, V, acf, Toep, prior = NULL, noSigma = FALSE,
                     calc.prior = TRUE, debug = FALSE) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- toep.suff(Y = Y, X = X, V = V, acf = acf, Toep = Toep)
  }
  n <- suff$n
  Betahat <- suff$Beta.hat
  T <- suff$T
  S <- suff$S
  noBeta <- is.null(T)
  p <- ifelse(noBeta, 0, nrow(Betahat))
  q <- nrow(S)
  # get prior
  if(is.null(prior) ||
     !all(sort(names(prior)) == .PriorNames) ||
     any(sapply(prior, is.null))) {
    prior <- .DefaultPrior(prior, p, q, noSigma)
  }
  Lambda <- prior$Lambda
  Omega <- prior$Omega
  Psi <- prior$Psi
  nu <- prior$nu
  # change noBeta and noSigma to make consistent with prior
  if(!noBeta) noBeta <- (length(Omega) == 1) && is.na(Omega)
  if(!noSigma) noSigma <- is.na(nu)
  if(noBeta && noSigma) {
    warning("Beta and Sigma both known.  calc.prior ignored.")
    return(list(Lambda = 0, Omega = NA, Psi = 0, nu = NA))
  }
  if(debug) browser()
  # calculate posterior parameters
  if(noBeta) {
    Lambdahat <- NA
    Omegahat <- NA
    noOmega <- TRUE
  } else {
    Omegahat <- T+Omega
    Lambdahat <- solve(Omegahat, T%*%Betahat + Omega%*%Lambda)
    noOmega <- all(Omega == 0)
  }
  if(noSigma) {
    Psihat <- S
    nuhat <- NA
  } else {
    Psihat <- Psi + S
    nuhat <- nu + (n-p*noOmega)
  }
  if(!noBeta) {
    Psihat <- Psihat + crossprod(Betahat, T%*%Betahat) +
      crossprod(Lambda,Omega%*%Lambda) -
        crossprod(Lambdahat,Omegahat%*%Lambdahat)
  }
  out <- list(Lambda = Lambdahat, Omega = Omegahat,
              Psi = Psihat, nu = nuhat)
  if(calc.prior) {
    prior <- list(Lambda = Lambda, Omega = Omega, Psi = Psi, nu = nu)
    out <- c(out, list(prior = prior))
  }
  out
}

if(FALSE) {
  Lambda <- prior$Lambda
  Omega <- prior$Omega
  nu <- prior$nu
  Psi <- prior$Psi
  if(is.null(nu)) nu <- 0
  if(is.null(Psi)) Psi <- 0
  if(!noBeta) {
    noBeta <- (length(Omega) == 1) && is.na(Omega)
  }
  if(!noSigma) {
    noSigma <- is.na(nu)
  }
  if(noBeta && noSigma) {
    warning("Beta and Sigma both known.  calc.prior ignored.")
    return(list(Lambda = 0, Omega = NA, Psi = 0, nu = NA))
  }
  if(noBeta) {
    Lambda <- 0
    Omega <- NA
    noOmega <- TRUE
  } else {
    if(is.null(Lambda)) Lambda <- matrix(0,p,q)
    Omega <- prior$Omega
    if(is.null(Omega)) Omega <- 0
    noOmega <- all(Omega == 0)
    if(noOmega) Omega <- matrix(0,p,p)
  }
  if(noSigma) {
    nu <- NA
  }
}
