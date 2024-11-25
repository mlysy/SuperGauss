#' Update the k-th factor process
#' @param k Integer for the factor index
#' @param Y `N x D` matrix of observations
#' @param X_minus_k `N x (K-1)` matrix of the current states of the factor processes
#' excluding the k-th factor process
#' @param sig2 Length `D` vector of diagonal elements of the noise process
#' @param beta `K x D` matrix of factor loadings
#' @param alpha Length `D` vector of bias parameters
#' @param theta Vector of factor process covariance kernel parameters
#' @param delta_t Length of each time bin
#' @param acf A function defining the acf of x
#' @param rand_sample Sample from the conditional density?
#' @param logdens Evaluate the log density of the full conditional? 
#' @return A vector of the sample of current x_k and/or the log density function
#' of the full conditional
#' @noRd
#' @export
cond_x_k <- function(k, Y, X_minus_k, sig2, beta, alpha, theta, delta_t, 
                     acf=rbf_acf, rand_sample=TRUE, logdens=FALSE){
  n_time <- dim(Y)[1]
  n_proc <- dim(Y)[2]
  # Compute acf for V_theta
  lag <- (0:(n_time-1))*delta_t
  V_acf <- do.call(function(...){acf(lag, ...)}, as.list(theta))
  
  Sig_inv_beta_k <- (1/sig2)*beta[k,]
  const <- sum(beta[k,] *Sig_inv_beta_k) # Scalar
  alpha_mat <- matrix(rep(alpha, n_time), nrow=n_time,ncol=n_proc, byrow = TRUE)
  yk_star <- Y - X_minus_k %*% beta[-k, , drop=F] - alpha_mat
  
  Q_acf <- const * V_acf
  Q_acf[1] <- Q_acf[1] + 1
  
  # Construct toeplitz objects
  #Q_tz <- toeplitz(Q_acf) # using the R implementation of topelitz class
  Q_tz <- Toeplitz$new(acf=Q_acf) # using the c++ implementation of toeplitz class
  V_tz <- Toeplitz$new(acf=V_acf)
  
  yk_Siginv_betak <- yk_star%*%Sig_inv_beta_k
  
  #------- log density of the full conditional----------
  logdens_function <- NULL
  if (logdens){
    VQinv <- V_tz$prod(Q_tz$solve(diag(N)))
    logdens_function <- function(x_k){
      mvtnorm::dmvnorm(x_k, mean=VQinv %*% yk_Siginv_betak,
                       sigma=VQinv, log=TRUE)
    }
  }
  
  #------- Generate a rv from the full conditional ---------
  out_sample <- NULL
  if (rand_sample){
    eps1 <- rnormtz(n=1, V_acf)
    eps2 <- rnorm(n=rep(1, n_time),
                  mean=yk_Siginv_betak, sd=rep(sqrt(const), n_time))
    #comp1 <- solve(Q_tz, eps1) # using the R implementation 
    comp1 <- Q_tz$solve(eps1) # using the c++ implementation
    #comp2 <- Vx_prod(V_acf, solve(Q_tz, eps2)) # using the R implementation
    comp2 <- V_tz$prod(Q_tz$solve(eps2)) # using the c++ implementation
    out_sample <- comp1 + comp2
  }
  return(list(sample=out_sample, logdens=logdens_function))
}

#' Alternative update noise variance matrix all together
#' @details Arguments are similar to those in `cond_sig2()`
#' @noRd
cond_sig2 <- function(Y, X, beta, alpha, sig_alpha, sig_beta,
                      rand_sample = TRUE, logdens = FALSE){
  n_time <- dim(Y)[1]
  n_proc <- dim(Y)[2]
  alpha_mat <- matrix(rep(alpha, n_time), nrow=n_time,ncol=n_proc, byrow = TRUE)
  y_minus <- Y - X%*%beta - alpha_mat
  if (length(sig_alpha)) sig_alpha <- rep(sig_alpha, n_proc)
  if (length(sig_beta)) sig_beta <- rep(sig_beta, n_proc)
  sig2_star <- rep(0, n_proc)
  gamma_shape <- function(d) sig_alpha[d] + .5*n_time
  gamma_scale <- function(d) sig_beta[d] + .5*crossprod(y_minus[,d])
  
  logdens_function <- NULL
  if (logdens){
    logdens_function <- function(sig2){
      sum(sapply(1:n_proc, 
                 function(d) {
                   dgamma(x=1/sig2[d], shape=gamma_shape(d), 
                          rate=gamma_scale(d), log=TRUE) - 2*log(sig2[d])
                   
                   }))
    }
  }
  out_sample <- NULL
  if (rand_sample){
    for (d in 1:n_proc){
      sig2_star[d] <- 1/rgamma(1, shape=gamma_shape(d),
                               rate=gamma_scale(d))
    }
    out_sample <- sig2_star
  }

  return(list(sample=out_sample, logdens=logdens_function))
}

#' Update of noise variance matrix element by element
#' @param y_d Length `N` vector of the `d`-th observed processes
#' @param X `N x K` (`n_time x n_factor`) matrix of factor processes
#' @param beta_d `K x 1` matrix of factor loadings
#' @param alpha_d Scalar `d`-th bias parameter
#' @param sig_alpha Scalar shape parameter for the prior on sigma_d
#' @param sig_beta Scalar scale parameter for the prior on sigma_d
#' @param rand_sample Sample from the conditional density?
#' @param logdens Evaluate the log density of the full conditional? 
#' @return Length `D` samples of sig2 (variances and/or the log density of the 
#' full conditional
#' @noRd
#' @export
cond_sig2_d <- function(y_d, X, beta_d, alpha_d, sig_alpha_d, sig_beta_d,
                        rand_sample = TRUE, logdens = FALSE){
  n_time <- length(y_d)
  y_minus <- y_d - X%*%beta_d - alpha_d
  gamma_shape <- sig_alpha_d + .5*n_time
  gamma_scale<- sig_beta_d + .5*crossprod(y_minus)
  
  logdens_function <- NULL
  if (logdens){
    logdens_function <- function(sig2_d){
      dgamma(x=1/sig2_d, shape=gamma_shape, 
             rate=gamma_scale, log=TRUE) - 2*log(sig2_d)
    }
  }
  
  out_sample <- NULL
  if (rand_sample) out_sample <- 1/rgamma(1, shape=gamma_shape, rate=gamma_scale)
  return(list(sample=out_sample, logdens=logdens_function))
}

#' Update the d-th column of factor loading matrix beta
#' @param y_d Length `N` vector of the `d`-th observed processes
#' @param X `N x K` (`n_time x n_factor`) matrix of factor processes
#' @param XtX `K x K` matrix of `crossprod(X)`
#' @param sig_d `d`-th diagonal element of observation cov matrix `Sigma`
#' @param alpha_d `d`-th element of the bias parameter `alpha`
#' @param psi_d Length `K` mean vector of the Normal prior on `beta_d`
#' @param S_d `K x K` covariance matrix of the Normal prior on `beta_d`
#' @param rand_sample Sample from the conditional density?
#' @param logdens Evaluate the log density of the full conditional? 
#' @return Vector sample of beta and/or the log density of the full conditional
#' @noRd
#' @export
cond_beta_d <- function(y_d, X, XtX, sig_d, alpha_d, psi_d, S_d,
                        rand_sample = TRUE, logdens = FALSE){
  K <- nrow(S_d)
  S_d_chol <- chol(S_d)
  S_d_inv <- backsolve(r = S_d_chol,
                       x = backsolve(r = S_d_chol, x = diag(K), transpose = T))
  S_d_inv_psi_d <- backsolve(r = S_d_chol,
                             x = backsolve(r = S_d_chol, x = psi_d, transpose = T))
  S_star <- solve(S_d_inv + XtX/sig_d)
  Psi_star <- S_star %*% (S_d_inv_psi_d + crossprod(X, y_d-alpha_d)/sig_d)
  
  logdens_function <- NULL
  if (logdens) {
    logdens_function <- function(beta_d){
      mvtnorm::dmvnorm(beta_d, mean=Psi_star, sigma=S_star, log=TRUE)
    }
  }
  
  out_sample <- NULL
  if (rand_sample) out_sample <- MASS::mvrnorm(n=1, mu=Psi_star, Sigma=S_star)
  
  return(list(sample=out_sample, logdens=logdens_function))
}

#' Update of bias parameter element by element
#' @param d Integer
#' @param y_d Length `N` vector of the `d`-th observed processes
#' @param X `N x K` (`n_time x n_factor`) matrix of factor processes
#' @param sig_d Scalar `d`-th diagonal element of observation cov matrix `Sigma`
#' @param beta_d `K x 1` matrix of factor loadings
#' @param nu_d Scalar mean parameter for the prior on alpha_d
#' @param tau_d Scalar variance parameter for the prior on alpha_d
#' @param rand_sample Sample from the conditional density?
#' @param logdens Evaluate the log density of the full conditional? 
#' @return Scalar sample of the `d`-th element of alpha and/or the log density  
#' of the full conditional
#' @noRd
#' @export
cond_alpha_d <- function(y_d, X, sig_d, beta_d, nu_d, tau_d,
                         rand_sample = TRUE, logdens = FALSE){
  n_time <- length(y_d)
  y_minus <- y_d - X%*%beta_d
  var_star <- 1/(1/tau_d + n_time/sig_d)
  sd_star <- sqrt(var_star)
  mean_star <- var_star * (nu_d/tau_d + sum(y_minus)/sig_d)
  
  logdens_function <- NULL
  if (logdens){
    logdens_function <- function(alpha_d){
      dnorm(alpha_d, mean=mean_star, sd=sd_star, log=TRUE)
    }
  }
  
  out_sample <- NULL
  if (rand_sample) out_sample <- rnorm(1, mean=mean_star, sd=sd_star)
  return(list(sample=out_sample, logdens=logdens_function))
}

#' Update the k-th theta using a Metropolis-Hastings step
#' @param x_k Length `N` vector of the `k`-th factor process.
#' @param current_theta Length 2 vector of current values of `current_theta`.
#' @param delta_t Length of each time bin.
#' @param acf The autocorrelation function for the `k`-th factor process.
#' @param ell_alpha Scalar shape parameter for the IG prior on theta2_k, where
#' theta2 is the lengthscale, if applicable.
#' @param ell_beta Scalar scale parameter for the IG prior on theta2_k.
#' If either of `ell_alpha` or `ell_beta` is NULL or NA, a uniform prior is used, 
#' but this is not recommended as a soft regularizing prior is need so that the
#' GP lengthscale does not shrink to zero.
#' @param prop_var_scale Vector of same length as current theta for the proposal
#' distribution SD scale. If is a scalar, the same scalar is used for proposing 
#' for all thetas.
#' @param fix_gp_widthscale Is the widthscale fixed for each factor process? 
#' If TRUE, the the widthscale values are fixed to be the first row of `theta`.
#' @param additional_prior Additional prior function to be applied on theta.
#' @return Length 2 list of `theta` and `accept`.
#' @noRd
#' @export
cond_theta_k <- function(x_k, current_theta, delta_t, acf=rbf_acf,
                         ell_alpha=NULL, ell_beta=NULL,
                         prop_var_scale=NULL,
                         fix_gp_widthscale=FALSE,
                         additional_prior=NULL){
  n_theta <- length(current_theta)
  if (is.null(prop_var_scale)) {
    prop_var_scale <- rep(.5*delta_t, length(current_theta))
  } else if (length(prop_var_scale)==1){
    prop_var_scale <- rep(prop_var_scale, n_theta)
  }
  n_time <- length(x_k)
  lag <- (0:(n_time-1))*delta_t
  current_acf <- do.call(function(...){acf(lag, ...)}, as.list(current_theta))
  current_lpdf <- dnormtz(X=x_k, acf=current_acf, log=T, method="gschur")
  if (is.na(current_lpdf)){
    stop("Cannot evaluate GP density due to ill-conditioned covariance matrix.")
  } 
  if (fix_gp_widthscale){
    new_theta <- c(current_theta[1],
                   rnorm(length(current_theta[-1]), current_theta[-1], prop_var_scale[-1]))
    lproposal <- function(x){
      sum(dnorm(x, current_theta[-1], prop_var_scale[-1], log=TRUE))
    }
  } else{
    new_theta <- rnorm(rep(1, n_theta), current_theta, prop_var_scale)
    lproposal <- function(x){
      sum(dnorm(x, current_theta, prop_var_scale, log=TRUE))
    }
  }
  if (is.null(ell_alpha) | is.null(ell_beta) | is.na(ell_alpha) | is.na(ell_beta)){
    new_prior <- 0
    current_prior <- 0
  } else{
    new_prior <- dgamma(x=1/new_theta[2], shape=ell_alpha, 
                        rate=ell_beta, log=TRUE) - 2*log(new_theta[2])
    current_prior <- dgamma(x=1/current_theta[2], shape=ell_alpha,
                            rate=ell_beta, log=TRUE) - 2*log(new_theta[2])
  }
  
  # Additional prior
  if (!is.null(additional_prior)){
    new_prior <- new_prior + additional_prior(new_theta)
    current_prior <- current_prior + additional_prior(current_theta)
  }
  
  new_acf <- do.call(function(...){acf(lag, ...)}, as.list(new_theta))
  new_lpdf <- dnormtz(X=x_k, acf=new_acf, log=T, method="gschur")
  if (is.na(new_lpdf)) new_lpdf <- -Inf
  if ((new_lpdf+new_prior-current_lpdf-current_prior) > log(runif(1))){
    list(current=new_theta, accept=1, lproposal=lproposal,
         new_lpdf=new_lpdf, current_lpdf=current_lpdf)
  } else{
    list(current=current_theta, accept=0, lproposal=lproposal,
         new_lpdf=new_lpdf, current_lpdf=current_lpdf)
  }
}


#' Gibbs sampler for the GPFA
#' @param n_iter Number of iterations.
#' @param Y `N x D` (`n_time x n_proc`) matrix of observed processes.
#' @param init_param A list of parameters:
#' - X: `N x K` (`n_time x n_factor`) matrix of factor processes.
#' - sig2: length `D` vector of diagonal elements (variances) of observation cov matrix `Sigma`.
#' - beta: `K x D` matrix of coefficients.
#' - alpha: Length `D` vector of bias parameters
#' - theta: A length `K` list of vectors of hyperparameters for the acf function for the
#' factor processes, where first element of vector is the widthscale (var), second element
#' of the vector is the lengthscale, and rest are other optional ones specific to the
#' `k`-th kernel function.
#' @param prior_list A list the following prior hyperparameters:
#' - Length `D` vectors `sig_alpha` and `sig_beta`: the shape and rate parameters
#' for the IGamma priors on `(sig2_1, ..., sig2_D)`. If a scalar is provided,
#' then the same shape/rate is used on all priors.
#' - Length `K` vector `psi_d`: the mean of the Normal prior on `beta_d`.
#' - `K x K` matrix `S_d`: the covariance matrix of the Normal prior on `beta_d`.
#' - Length `D` vector `nu`: the mean of the Normal prior on each element of `alpha`.
#' - Length `D` vector `tau`: the sd of the Normal prior on each element of `alpha`.
#' - Length `K` vectors `ell_alpha` and `ell_beta`: the shape and rate parameters 
#' for the IGamma priors on `(lengthscale_1, ..., lengthscale_K)` of the GP 
#' parameters. If a scalar is provided, then the same shape/rate is used on all 
#' priors.
#' @param delta_t Length of each time bin.
#' @param acf_list A list of autocorrelation functions for the `K` factor processes.
#' Default is the ACF of the radial basis kernel.
#' @param theta_prop_scale List of vectors: scales of the variance for proposal 
#' for each element of theta.
#' @param fix_gp_widthscale Is the widthscale fixed for each factor process? 
#' If TRUE, the the widthscale values are fixed to be the first row of `theta`.
#' @param print_every Integer. Print progress and acceptance rate every `print_every`
#' steps.
#' @param fix_param Optionally, parameters can be fixed in sampling. This
#' argument is a named list where 
#' - `names(fix_param)` is a subset of `names(init_param)`,
#' - `length(fix_param$p)` equals `length(init_param$p)` for a parameter "p",
#' - "p" entries with NAs are fixed in sampling.
#' This argument is only used for development purposes for now and is not 
#' fully tested. When fixing parameter values, can only fix all values of a column
#' of `X` and `beta`, or all theta values for a latent path 
#' (i.e., all of `theta[[k]]`).
#' @return Samples of parameters.
#' @export
gpfa_gibbs_sampler <- function(n_iter, Y, init_param, prior_list, delta_t, 
                               acf_list=rep(list(rbf_acf), nrow(init_param$X)),
                               theta_prop_scale=NULL, fix_gp_widthscale=FALSE, 
                               print_every=10, fix_param=list()){
  n_time <- dim(Y)[1]
  n_proc <- dim(Y)[2]
  X <- init_param$X
  sig2 <- init_param$sig2
  beta <- init_param$beta
  alpha <- init_param$alpha
  theta <- init_param$theta
  n_factor <- dim(X)[2]
  X_sam <- array(NA, dim=c(n_iter, dim(X)))
  sig2_sam <- matrix(NA, nrow=n_iter, ncol=n_proc)
  beta_sam <- array(NA, dim=c(n_iter, dim(beta)))
  alpha_sam <- matrix(NA, nrow=n_iter, ncol=n_proc)
  theta_sam <- array(NA, dim=c(n_iter, length(unlist(theta))))
  thetak_accept <- rep(0, n_factor)
  for (i in 1:n_iter){
    #---------- Update X -----------------
    for (k in 1:n_factor){
      if (is.null(fix_param$X) | !all(is.na(fix_param$X[, k]))){
        X[, k] <- cond_x_k(k, Y, X[,-k, drop=F], sig2, beta, alpha, theta[[k]], 
                           delta_t, acf=acf_list[[k]])$sample
      }
    }
    X_sam[i,,] <- X
    #---------- Update variance -------------
    for (d in 1:n_proc){
      if (is.null(fix_param$sig2) | !all(is.na(fix_param$sig2[d]))){
        sig2[d] <- cond_sig2_d(Y[,d], X, beta[,d], alpha[d],
                               prior_list$sig_alpha[d], 
                               prior_list$sig_beta[d])$sample
      }
    }
    sig2_sam[i,] <- sig2
    #---------- Update beta --------------
    XtX <- crossprod(X)
    for (d in 1:n_proc){
      if (is.null(fix_param$beta) | !all(is.na(fix_param$beta[,d]))){
        beta[,d] <- cond_beta_d(Y[,d], X, XtX, sig2[d], alpha[d],
                                prior_list$psi_d, prior_list$S_d)$sample
      }
    }
    beta_sam[i,,] <- beta
    #---------- Update alpha --------------
    for (d in 1:n_proc){
      if (is.null(fix_param$alpha) | !all(is.na(fix_param$alpha[d]))){
        alpha[d] <- cond_alpha_d(Y[,d], X, sig2[d], beta[,d], 
                                 prior_list$nu[d], prior_list$tau[d])$sample
      }
    }
    alpha_sam[i,] <- alpha
    #---------- Update theta -----------------
    for (k in 1:n_factor){
      if (is.null(fix_param$theta) | !all(is.na(fix_param$theta[[k]]))){
        theta_update <- cond_theta_k(X[,k], theta[[k]], delta_t, acf_list[[k]],
                                     prior_list$ell_alpha[k], 
                                     prior_list$ell_beta[k],
                                     theta_prop_scale[[k]], fix_gp_widthscale,
                                     prior_list$theta_additional_prior)
        theta[[k]] <- theta_update$current
        thetak_accept[k] <- thetak_accept[k] + theta_update$accept
      }
    }
    theta_sam[i,] <- unlist(theta)
    if((print_every > 0) && (i %% print_every == 0)) {
      cat("Iteration", i, "-- theta acc_rate:", signif(thetak_accept/i, 2), "\n")
    }
  }
  list(X_sam=X_sam, sig2_sam=sig2_sam, beta_sam=beta_sam, alpha_sam=alpha_sam,
       theta_sam=theta_sam, theta_accept=thetak_accept/n_iter)
}
