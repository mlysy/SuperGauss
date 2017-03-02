#' Calculate the sufficient statistics of an LMN model using Superfast-Toeplitz.
#'
#' @param Y an (n x q) matrix.
#' @param X an (n x p) matrix. can also be passed as
#' \itemize{
#'    \item \code{X = 0}: in which case there is no intercept
#'    \item \code{X != 0}: in which case a scaled intercept X = X * matrix(1, n, 1)
#'    is assumed.
#' }
#' @param V Either: (1) an \code{n x n} full matrix, (2) an vector of length \code{n} such that \code{V = diag(V)}, (3) a scalar, such that \code{V = V * diag(n)}.
#' @param acf a vector of length n such that V = toeplitz(acf).
#' @param npred: an integer.  if nonzero, returns sufficient statistics to make predictions for \code{y | Y, Beta, Sigma} (see details).
#' @details
#' The model is defined as
#' \deqn{Y ~ MNorm(X*Beta, V, Sigma),}
#' where MNorm is the Matrix-Normal distribution, i.e.
#' \deqn{vec(Y) ~ N( vec(X*Beta), Sigma \otimes V ).}
#'
#' The log-likelihood is
#' \deqn{ll(Beta, Sigma) = (-1/2) * ( trace( Sigma^{-1} ( S + t(Beta-B) T (Beta-B)) ) + n log(|Sigma|) + q ldV ).}
#' When \code{T = 0} and \code{S = 0} respectively, \code{Beta = 0} or \code{Sigma = diag(q)} is assumed to be known.
#'
#' If \code{npred > 0}, summary statistics are produced for the prediction of \code{y} in
#'    \deqn{\code{YY = rbind(Y, y) ~ MNorm(XX*Beta, VV, Sigma)},}
#'    where \code{VV = rbind(cbind(V, w), cbind(t(w), v))} and \code{XX = rbind(X, x)}.
#'    In this case the inputs are \code{Y}, \code{V = VV} and \code{X = XX}. The predictive distribution for new observations y is
#' \deqn{y | Y ~ MNorm(Ap + Xp * Beta, Vp, Sigma).}
#' \code{Xp = 0} indicates that \code{Beta = 0} is known.
#'
#' @return a list with the following elements:
#' \itemize{
#'    \item \code{Betahat = (t(X)V^{-1}X)^{-1}t(X)V^{-1}Y}
#'    \item \code{T = t(X)V^{-1}X}
#'    \item \code{S = t(Y-X*Betahat)V^{-1}(Y-X*Betahat)}
#'    \item \code{ldV = log(|V|)}
#'    \item \code{n = nrow(Y)}
#'    \item \code{noBeta} and \code{noSigma}
#' }
#' In addition, when \code{npred > 0}:
#' \itemize{
#'    \item \code{Ap = t(w)V^{-1}Y}
#'    \item \code{Xp = x - t(w)V^{-1}X}
#'    \item \code{Vp = v - t(w)V^{-1}w}
#' }
#' @export
toep.suff <- function(Y, X, V, acf, npred = 0,
                     Toep, debug = FALSE) {
  # size of problem
  n <- nrow(Y)
  q <- ncol(Y)
  if(length(X) == 1) {
    X <- matrix(X, n+npred, 1)
  }
  noBeta <- all(X == 0)
  p <- (!noBeta) * ncol(X)
  # dimension check for Toep
  if(missing(Toep)){
    Toep <- new(Toeplitz, n)
  }
  else{
    # assuming that we have a given environment Toep
    if(!(n == Toep$DimCheck())){
      stop("Given Toep is incompatible with nrow(Y).")
    }
  }
  # inner product calculations
  Z <- matrix(0, n, q+p+npred)
  Z[,1:q] <- Y
  if(!noBeta) {
    Z[,q+(1:p)] <- X[1:n,]
  }
  # variance type
  if(!missing(V)) {
    if((length(V) == 1) || (length(V) == n + npred)) {
      var.type <- "diag"
    } else if(is.matrix(V) && all(dim(V) == n + npred)) {
      var.type <- "full"
    } else {
      stop("Given V is incompatible with nrow(Y) and npred.")
    }
  } else if(!missing(acf)) {
    if(length(acf) == n + npred) {
      var.type <- "acf"
    } else {
      stop("Given acf is incompatible with nrow(Y) and npred.")
    }
  } else {
    stop("Either V or acf must be provided.")
  }
  if(debug) browser()
  if(var.type == "full") {
    C <- chol(V[1:n,1:n])
    if(npred > 0) {
      Z[,q+p+(1:npred)] <- V[1:n,n+(1:npred)]
    }
    # IP' * V^{-1} * IP
    IP <- crossprod(backsolve(r = C, x = Z, transpose = TRUE))
    # log|V|
    ldV <- 2 * sum(log(diag(C)))
  } else if(var.type ==  "diag") {
    IP <- matrix(0, q+p+npred, q+p+npred)
    if(length(V) == 1) {
      IP[1:(q+p),1:(q+p)] <- crossprod(Z[,1:(q+p)])/V
      ldV <- n * log(V)
      IP[1:(q+p),1:(q+p)] <- crossprod(Z[,1:(q+p)], Z[,1:(q+p)]/V[1:n])
      ldV <- sum(log(V[1:n]))
    }
  } else if(var.type == "acf") {
    if(npred > 0) {
      Z[,q+p+(1:npred)] <- toeplitz2(rev(acf[1+1:n]), acf[n+1:npred])
    }
    # replace this by GSchur IP-ldV calculation
    # Z' * toeplitz(acf[1:n])^{-1} * Z
    # if is.numeric(acf), create space
    # if class(acf) == "Toeplitz", don't create.
    Toep$AcfInput(acf[1:n])
    IP <- crossprod(Z, Toep$Solve(Z))
    ldV <- Toep$Det()
    # DL <- DurbinLevinsonEigen(X = t(Z), Y = matrix(0),
    #                           acf = acf[1:n], calcMode = 1)
    # IP <- DL$IP
    # ldV <- DL$ldV
    
  } else {
    stop("Unrecognized variance type.")
  }
  # convert inner products to sufficient statistics
  S <- IP[1:q,1:q,drop=FALSE]
  if(noBeta) {
    Beta.hat <- NULL
    T <- NULL
  } else {
    T <- IP[q+(1:p),q+(1:p),drop=FALSE]
    Beta.hat <- solve(T, IP[q+(1:p),1:q,drop=FALSE])
    S <- S - IP[1:q,q+(1:p)] %*% Beta.hat
  }
  # predictive distribution
  if(npred > 0) {
    Ap <- IP[q+p+(1:npred),1:q,drop=FALSE]
    if(noBeta) {
      Xp <- NULL
    } else {
      Xp <- X[n+(1:npred),] - IP[q+p+(1:npred),q+(1:p),drop=FALSE]
    }
    if(var.type == "full") {
      Vp <- V[n+(1:npred),n+(1:npred),drop=FALSE]
    } else if(var.type == "diag") {
      if(length(V) == 1) {
        Vp <- diag(V, npred)
      } else {
        Vp <- diag(V[n+(1:npred)], npred)
      }
    } else if(var.type == "acf") {
      Vp <- toeplitz(acf[1:npred])
    } else {
      stop("Unrecognized variance type.")
    }
    Vp <- Vp - IP[q+p+(1:npred),q+p+(1:npred)]
  }
  # output
  ans <- list(Beta.hat = Beta.hat, S = S, T = T, ldV = ldV, n = n)
  if(npred > 0) {
    ans <- c(ans, list(Ap = Ap, Xp = Xp, Vp = Vp))
  }
  ans
}
