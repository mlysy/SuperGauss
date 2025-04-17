require(MASS)
require(SuperGauss)
require(tidyverse)
require(ggpubr)
require(testthat)
sampler_scripts <- list.files("../", pattern="*.R$",
                              full.names=TRUE, ignore.case=TRUE)
sapply(sampler_scripts, source, .GlobalEnv)
# Notation:
#   y_n: 1 x n_proc observed spike trains in time bin n (n_time)
#   beta: n_factor x n_proc design matrix
#   x(t): 1 x n_factor process
#   Sigma: n_proc x n_proc diagonal matrix
#   x_k: 1 x n_time factor process binned realization for factor k
#   V^k: n_time x n_time temporal covariance matrix for x_k
# Model:
#   Vector: y_n | x(t) ~ N(x(n*delta_t)*beta + alpha, Sigma)
#   Vector: x_k ~ N(0, V^k)

#-------------- Simulate toy data for unit tests ----------------
set.seed(123)
N <- 2000
D <- 4
K <- 2
theta <- cbind(c(1, 0.5), c(1.2, 0.6))^2
beta <- matrix(rnorm(K*D), nrow=K, ncol=D)
delta_t <- 0.5
sig2 <- rnorm(D, 0.5, 0.5)^2
X <- sapply(1:K,
            function(k){
              x_acf <- rbf_acf((0:(N-1))*delta_t,
                               var=theta[1,k], ell=theta[2,k])
              rnormtz(n = 1, x_acf)
            })
alpha <- rnorm(D, mean=4, sd=1)
Y <- X %*% beta + matrix(rep(alpha, N), nrow=N, ncol=D, byrow = TRUE) +
  mvrnorm(n=N, mu=rep(0, D), Sigma=diag(sig2))
prior_list <- list(sig_alpha=rep(2, D),
                   sig_beta=rep(2, D),
                   psi_d=rep(1, K),
                   S_d=diag(K),
                   nu=rep(2, D),
                   tau=rep(25, D),
                   ell_alpha=rep(5,K),
                   ell_beta=rep(5,K))

#-------------- Test via Bayes rule --------
# Idea: p(theta' | data) / p(theta | data) = p(theta', data) / p(theta, data)
# Ref: Grosse and Duvenaud (2014) https://arxiv.org/pdf/1412.5218

source("helper_gibbs_unittests.R") # load all test functions

current_state <- list(X=X, sig2=sig2, beta=beta, alpha=alpha, theta=theta)
# Test conditional update of x_k
for (k in 1:K){
  test_cond_x_k(k, Y, current_state, prior_list, delta_t)
}

# Test conditional update of beta_d
for (d in 1:D){
  test_cond_beta_d(d, Y, current_state, prior_list, delta_t)
}

# Test conditional update of sig2_d
for (d in 1:D){
  test_cond_sig2_d(d, Y, current_state, prior_list, delta_t)
}

# Test conditional update of alpha_d
for (d in 1:D){
  test_cond_alpha_d(d, Y, current_state, prior_list, delta_t)
}

# Test MH block for updating theta_k
n_sim <- 1000
n_warmup <- 500
theta_sam_list <- vector("list", 2)
for (k in 1:K){
  theta_sam_list[[k]] <- test_cond_theta_k(k, true_theta_k=theta[,k],
                                  current_theta_k=c(theta[1,k],0.2), X_k=X[,k],
                                  delta_t=delta_t, prior_list,
                                  prop_var_scale=rep(0.02, 2),
                                  fix_gp_widthscale=TRUE,
                                  n_sim=n_sim, n_warmup=n_warmup)
}
