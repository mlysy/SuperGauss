#' Calculates the joint log likelihood at the current state of parameters.
#' @param Y Matrix of data.
#' @param current_state A named list of `X`, `sig2`, `beta`, `alpha`, `theta`
#' containing the current states of the model parameters.
#' @param prior_list A named list of `sig_alpha`, `sig_beta`, `psi_d`, `S_d`,
#' `nu`, `tau`, `ell_alpha`, `ell_beta`.
#' containing the prior values for the inverse gamma prior on `sig2` and the
#' Normal prior on `beta`.
#' @param delta_t Scalar value of bin size for the processes.
#' @return Scalar value of the joint log likelihood evaluated at the current state.
joint_log_prob <- function(Y, current_state, prior_list, delta_t){
  X <- current_state$X
  sig2 <- current_state$sig2
  beta <- current_state$beta
  alpha <- current_state$alpha
  theta <- current_state$theta

  logpdf <- 0
  N <- dim(Y)[1]
  D <- dim(Y)[2]
  K <- dim(X)[2]
  # Likelihood
  for (n in 1:N){
    logpdf <- logpdf + mvtnorm::dmvnorm(x=Y[n,],
                                        mean=as.vector(X[n,,drop=FALSE]%*%beta)+alpha,
                                        sigma=diag(sig2), log=TRUE)
  }

  # Factor processes
  lag <- (0:(N-1))*delta_t
  for (k in 1:K){
    acf_k <- rbf_acf(lag, theta[1,k], theta[2,k])
    logpdf <- logpdf + dnormtz(X=X[,k], acf=acf_k, log=TRUE, method="gschur")
  }

  # Priors
  for (d in 1:D){
    # beta prior
    logpdf <- logpdf + mvtnorm::dmvnorm(x=beta[,d], mean=prior_list$psi_d,
                                        sigma=prior_list$S_d, log=TRUE)
    # alpha prior
    logpdf <- logpdf + dnorm(x=alpha[d], mean=prior_list$nu[d],
                             sd=sqrt(prior_list$tau[d]), log=TRUE)

    # sigma2 prior
    logpdf <- logpdf +
      dgamma(x=1/sig2[d], shape=prior_list$sig_alpha[d],
             rate=prior_list$sig_beta[d], log=TRUE) - 2*log(sig2[d]) # invgamma
  }
  if (!is.null(prior_list$ell_alpha) & !is.null(prior_list$ell_beta)){
    for (k in 1:K){
      # theta prior
      logpdf <- logpdf +
        dgamma(x=1/theta[2,k], shape=prior_list$ell_alpha[k],
               rate=prior_list$ell_beta[k], log=TRUE) - 2*log(theta[2,k])
    }
  }
  logpdf
}

#' Test if the full conditional of x_k is doing what it is supposed to
#' @details Arguments are similar to those of `joint_log_prob()`.
test_cond_x_k <- function(k, Y, current_state, prior_list, delta_t){
  X <- current_state$X
  sig2 <- current_state$sig2
  beta <- current_state$beta
  alpha <- current_state$alpha
  theta <- current_state$theta

  # Current and new states
  xk_current <- X[, k]
  xk_new <- cond_x_k(k, Y, X[,-k, drop=F], sig2, beta, alpha, theta[,k],
                     delta_t, logdens=TRUE)
  X_new <- X
  X_new[,k] <- xk_new$sample
  new_state <- current_state
  new_state$X <- X_new

  # Full conditional evaluated at current and new states
  cond_xk <- xk_new$logdens
  cond_xk_current <- cond_xk(as.numeric(xk_current))
  cond_xk_new <- cond_xk(as.numeric(xk_new$sample))

  # Joint log likelihood at current and new states
  joint_lpdf_current <- joint_log_prob(Y, current_state, prior_list, delta_t)
  joint_lpdf_new <- joint_log_prob(Y, new_state, prior_list, delta_t)
  expect_equal(cond_xk_new-cond_xk_current, joint_lpdf_new-joint_lpdf_current)
}

#' Test if the full conditional of beta_d is doing what it is supposed to
#' @details Arguments are similar to those of `joint_log_prob()`.
test_cond_beta_d <- function(d, Y, current_state, prior_list, delta_t){
  X <- current_state$X
  sig2 <- current_state$sig2
  beta <- current_state$beta
  alpha <- current_state$alpha
  theta <- current_state$theta
  y_d <- Y[,d]
  XtX <- crossprod(X)
  # Current and new states
  bd_current <- beta[,d]
  bd_new <- cond_beta_d(y_d, X, XtX, sig2[d], alpha[d],
                        prior_list$psi_d, prior_list$S_d, logdens=TRUE)
  beta_new <- beta
  beta_new[,d] <- bd_new$sample
  new_state <- current_state
  new_state$beta <- beta_new

  # Full conditional evaluated at current and new states
  cond_bd <- bd_new$logdens
  cond_bd_current <- cond_bd(as.numeric(bd_current))
  cond_bd_new <- cond_bd(as.numeric(bd_new$sample))

  # Joint log likelihood at current and new states
  joint_lpdf_current <- joint_log_prob(Y, current_state, prior_list, delta_t)
  joint_lpdf_new <- joint_log_prob(Y, new_state, prior_list, delta_t)
  expect_equal(cond_bd_new-cond_bd_current, joint_lpdf_new-joint_lpdf_current)
}

#' Test if the full conditional of beta_d is doing what it is supposed to
#' @details Arguments are similar to those of `joint_log_prob()`.
test_cond_alpha_d <- function(d, Y, current_state, prior_list, delta_t){
  X <- current_state$X
  sig2 <- current_state$sig2
  beta <- current_state$beta
  theta <- current_state$theta
  y_d <- Y[,d]
  XtX <- crossprod(X)
  # Current and new states
  ad_current <- alpha[d]
  ad_new <- cond_alpha_d(y_d, X, sig2[d], beta[,d],
                        prior_list$nu[d], prior_list$tau[d], logdens=TRUE)
  alpha_new <- alpha
  alpha_new[d] <- ad_new$sample
  new_state <- current_state
  new_state$alpha <- alpha_new

  # Full conditional evaluated at current and new states
  cond_ad <- ad_new$logdens
  cond_ad_current <- cond_ad(as.numeric(ad_current))
  cond_ad_new <- cond_ad(as.numeric(ad_new$sample))

  # Joint log likelihood at current and new states
  joint_lpdf_current <- joint_log_prob(Y, current_state, prior_list, delta_t)
  joint_lpdf_new <- joint_log_prob(Y, new_state, prior_list, delta_t)
  expect_equal(cond_ad_new-cond_ad_current, joint_lpdf_new-joint_lpdf_current)
}

#' Test if the full conditional of sig2_d is doing what it is supposed to
#' @details Arguments are similar to those of `joint_log_prob()`.
test_cond_sig2_d <- function(d, Y, current_state, prior_list, delta_t){
  X <- current_state$X
  sig2 <- current_state$sig2
  beta <- current_state$beta
  alpha <- current_state$alpha
  theta <- current_state$theta
  # Current and new states
  sig2_d_new <- cond_sig2_d(Y[,d], X, beta[,d], alpha[d], prior_list$sig_alpha[d],
                            prior_list$sig_beta[d], logdens=TRUE)
  sig2_new <- sig2
  sig2_new[d] <- sig2_d_new$sample
  new_state <- current_state
  new_state$sig2 <- sig2_new

  # Full conditional evaluated at current and new states
  cond_sig2 <- sig2_d_new$logdens
  cond_sig2_current <- cond_sig2(as.numeric(sig2[d]))
  cond_sig2_new <- cond_sig2(as.numeric(sig2_d_new$sample))

  # Joint log likelihood at current and new states
  joint_lpdf_current <- joint_log_prob(Y, current_state, prior_list, delta_t)
  joint_lpdf_new <- joint_log_prob(Y, new_state, prior_list, delta_t)
  expect_equal(cond_sig2_new-cond_sig2_current, joint_lpdf_new-joint_lpdf_current)
}

#' Test if the Metropolis-Hastings update block is doing what it is supposed to
#' @param k Integer index
#' @param true_theta_k Length 2 vector
#' @param current_theta_k Length 2 vector
#' @param X_k Length `N` vector
#' @param prior_list A list of priors parameters.
#' @param delta_t Scalar value of bin size for the processes.
#' @param prop_var_scale Vector of same length as current theta for the proposal
#' distribution SD scale.
#' @param fix_gp_widthscale Is the widthscale fixed for each factor process?
#' If TRUE, the the widthscale values are fixed to be the first row of `theta`.
#' @param n_sim Number of iterations for M-H sampling`X`
#' @param print_every Integer for printing sampling progress
test_cond_theta_k <- function(k, true_theta_k, current_theta_k, X_k, delta_t,
                              prior_list,
                              prop_var_scale=delta_t/2, fix_gp_widthscale=FALSE,
                              n_sim=1e4, n_warmup=floor(n_sim/2),
                              print_every=n_sim/10){
  thetak_sam <- matrix(0, nrow=n_sim, ncol=2)
  n_acc <- 0
  for (iter in 1:n_sim){
    thetak_update <- cond_theta_k(X_k, current_theta_k,
                                  delta_t, ell_alpha=prior_list$ell_alpha[k],
                                  ell_beta=prior_list$ell_beta[k],
                                  prop_var_scale=prop_var_scale,
                                  fix_gp_widthscale=fix_gp_widthscale)
    thetak_sam[iter,] <- thetak_update$current
    current_theta_k <- thetak_update$current
    n_acc <- n_acc + thetak_update$accept
    if (iter %% print_every == 0) {
      cat("Iteration", iter, "-- acceptance rate:", n_acc/iter, "\n")
    }
  }
  thetak_sam <- thetak_sam[-(1:n_warmup),]
  par(mfrow=c(1,2))
  for (ii in 1:2){
    plot(density(thetak_sam[,ii]), xlab="theta value",
         main=paste0("theta",k,"_",ii))
    abline(v=true_theta_k[ii], col="red", cex=1.2)
  }
  return(thetak_sam)
}
