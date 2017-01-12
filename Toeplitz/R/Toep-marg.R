#' Marginal Log-Posterior for the LMN Model using Superfast Algorithm
#'
#' @param suff The sufficient statistics to be combined with the MNIW conjugate prior.
#' @param post The parameters of the conditional MNIW distribution.  If missing will use \code{prior} and \code{noSigma} to calculate.
#' @param prior The parameters of the prior.  These are required for correct normalization.
#' @param noSigma Used to calculate \code{post} if it is missing.
#' @export
toep.marg <- function(suff, Y, X, V, acf, Toep, post, prior, noSigma,
                     debug = FALSE) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- toep.suff(Y = Y, X = X, V = V, acf = acf, Toep = Toep)
  }
  n <- suff$n
  Betahat <- suff$Beta.hat
  S <- suff$S
  ldV <- suff$ldV
  noBeta <- is.null(suff$T)
  p <- ifelse(noBeta, 0, nrow(Betahat))
  q <- nrow(S)
  # calculate prior and posterior
  if(debug) browser()
  if(missing(prior)) {
    if(missing(post)) {
      post <- toep.post(suff = suff, noSigma = noSigma, Toep = Toep, prior = NULL)
      prior <- post$prior
    } else {
      prior <- post$prior
      if(is.null(prior)) {
        stop("post supplied with unspecified prior.")
      }
    }
  }
  if(is.null(prior) ||
     !all(sort(names(prior)) == sort(.PriorNames)) ||
     any(sapply(prior, is.null))) {
    prior <- .DefaultPrior(prior, p, q, noSigma)
  } else {
    noSigma <- is.na(prior$nu)
  }
  if(missing(post)) {
    post <- toep.post(suff = suff, prior = prior, noSigma = noSigma, Toep = Toep,
                     calc.prior = FALSE)
  }
  # posterior MNIW parameters
  Omegahat <- post$Omega
  Psihat <- post$Psi
  nuhat <- post$nu
  noBeta <- (length(Omegahat) == 1) && is.na(Omegahat)
  noSigma <- is.na(nuhat)
  # prior MNIW parameters
  nu <- prior$nu
  Psi <- prior$Psi
  Omega <- prior$Omega
  if(noBeta) {
    noOmega <- TRUE
  } else {
    noOmega <- all(Omega == 0)
  }
  if(noSigma) {
    noPsi <- TRUE
  } else {
    noPsi <- all(Psi == 0)
  }
  # marginal log-posterior
  lp <- 0
  if(!noBeta) {
    lp <- lp + ldet(Omegahat)
  }
  lp <- q * (lp + ldV)
  if(!noSigma) {
    lp <- -.5 * (nuhat * (ldet(Psihat) - q*log(2)) + lp)
    lp <- lp + lmgamma(.5*nuhat, q)
  } else {
    lp <- -.5 * (lp + sum(diag(Psihat)))
  }
  # normalize by prior
  lpi <- 0
  if(!noBeta && !noOmega) {
    lpi <- lpi + ldet(Omega)
  }
  lpi <- q * (lpi - (n - p*noOmega) * log(2*pi))
  if(!noSigma && !noPsi) {
    lpi <- -.5 * (nu * (ldet(Psi) - q*log(2)) + lpi)
    lpi <- lpi + lmgamma(.5*nu, q)
  } else {
    lpi <- -.5 * lpi
  }
  lp - lpi
}
