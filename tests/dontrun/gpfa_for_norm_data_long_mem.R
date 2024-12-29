require(MASS)
require(SuperGauss)
require(doParallel)
require(abind)
require(tidyverse)
require(ggpubr)
sampler_scripts <- list.files("gibbs_sampler", pattern="*.R$", 
                              full.names=TRUE, ignore.case=TRUE)
sapply(sampler_scripts, source, .GlobalEnv)

# Set what to do
do_sampling <- FALSE
do_sample_processing <- TRUE

# Set what to save
save_figs <- TRUE
save_samples_as_rds <- TRUE

# Set MCMC parameters
n_chain <- 4
n_sam <- 2e4
n_sam_batch_size <- n_sam/5

# Set simulation parameters
fix_theta1 <- TRUE
set.seed(1234)
n_proc <- 30
n_factor <- 3
n_time <- 2e4
delta_t <- 0.1
Sigma_diag <- rnorm(n_proc, 0.5, 0.1)^2
proc_per_group <- n_proc/n_factor
beta <- rbind(c(runif(proc_per_group, 3, 5), runif(proc_per_group, 0, 1), runif(proc_per_group, 0, 1)),
              c(runif(proc_per_group, 0, 1), runif(proc_per_group, 3, 5), runif(proc_per_group, 0, 1)),
              c(runif(proc_per_group, 0, 1), runif(proc_per_group, 0, 1), runif(proc_per_group, 3, 5)))
alpha <- rnorm(n_proc, mean=7, sd=2)
theta <- list(c(0.05, 0.05),
              c(0.05, 0.1),
              c(0.05, 0.12, 0.4))

#------- Simulate factor processes ----------
# Short memory latent process
acf1 <- rbf_acf((0:(n_time-1))*delta_t, var=theta[[1]][1], ell=theta[[1]][2])
acf2 <- rbf_acf((0:(n_time-1))*delta_t, var=theta[[2]][1], ell=theta[[2]][2])
# Long memory latent process
acf3 <- rquad_acf((0:(n_time-1))*delta_t, var=theta[[3]][1], k=theta[[3]][2], d=theta[[3]][3])
n_lags_to_plot <- 2e4

# # Plot ACFs
# acf_df <- data.frame(acf=c(acf1[1:n_lags_to_plot], acf2[1:n_lags_to_plot], acf3[1:n_lags_to_plot]), 
#                      lag=rep(0:(n_lags_to_plot-1), 3),
#                      type=c(rep("Short memory 1", n_lags_to_plot),
#                             rep("Short memory 2", n_lags_to_plot),
#                             rep("Long memory", n_lags_to_plot)))
# ggplot(data=acf_df)+
#   geom_line(aes(x=lag, y=acf))+facet_grid(rows=vars(type))+
#   labs(x="Lag", y="ACF")
# if (save_figs) ggsave("ggpfa-acf.pdf", width=10, height=6)

x1 <- rnormtz(n=1, acf1)
x2 <- rnormtz(n=1, acf2)
x3 <- rnormtz(n=1, acf3)
x_sim <- cbind(x1, x2, x3)
# # Plot x
# xk_titles <- factor(c("Short memory X1", "Short memory X2", "Long memory X3"),
#                     levels=c("Short memory X1", "Short memory X2", "Long memory X3"))
# x_plot_df <- data.frame(latent_x=c(x1, x2, x3), 
#                         time=rep(1:n_time, times=n_factor),
#                         type=rep(xk_titles, each=n_time))
# ggplot(data=x_plot_df)+
#   geom_line(aes(x=time, y=latent_x))+facet_grid(rows=vars(type))+
#   labs(x="Time", y="Simulated latent path")
# if (save_figs) ggsave("ggpfa-simulated-x.pdf", width=10, height=6)

# Simulate observed data
y_sim <- x_sim %*% beta + 
  matrix(rep(alpha, n_time), nrow=n_time, ncol=n_proc, byrow = TRUE) +
  mvrnorm(n=n_time, mu=rep(0, n_proc), Sigma=diag(Sigma_diag))
cat("Number of simulated Y_i smaller than 0 is:\n",
    sum(y_sim < 0), "\n")
# Define priors and random initial values
# log pi(d) = log (d + .5) + log(.5 - d)
prior_list <- list(sig_alpha=rep(1, n_proc),  # length D
                   sig_beta=rep(1, n_proc),   # length D
                   psi_d=rep(1, n_factor),      # length K
                   S_d=diag(rep(25,n_factor)),
                   nu=rep(0,n_proc),
                   tau=rep(100,n_proc),
                   ell_alpha=rep(NA, 3),
                   ell_beta=rep(NA, 3),
                   theta_additional_prior=list(NULL, NULL,
                     function(theta){
                       # Only apply this prior on the rf kernel d parameter
                       log(theta[3]+0.5)+log(0.5-theta[3])
                       }
                   ))

#-------- Initialize parameters for MCMC ---------------
# Good initial values: true values
init_param <- list(X=x_sim,
                   beta=beta,
                   alpha=alpha,
                   sig2=Sigma_diag,
                   theta=theta)

# # Very bad initial values: random
# init_x <- matrix(rep(0, n_time*n_factor), nrow=n_time, ncol=n_factor)
# for (i in 1:n_factor){
#   init_x[,i] <- ksmooth(1:n_time, y_sim[,i], "normal", bandwidth = 20)$y-alpha[i]
# }
# init_param <- list(X=init_x,
#                    beta=matrix(1, nrow=n_factor, ncol=n_proc),
#                    alpha=rep(1, n_proc),
#                    sig2=rep(1, n_proc),
#                    theta=list(rep(0.05, 2), rep(0.05, 2), rep(0.05, 3)))

# Initialize using factanal()
# init_param <- gpfa_init_params(y_sim, n_factor=n_factor, 
#                                acf_list=list(rbf_acf, rbf_acf, rquad_acf))

fix_param <- list() # any parameters to fix during sampling can go here
# Tune theta_prop_scale
theta_prop_scale <- list(c(0, 0.0008),  # first element set to 0 because fixing widthscale
                         c(0, 0.0008),
                         c(0, 0.001, 0.008))
# test <- gpfa_gibbs_sampler(100, y_sim, init_param, prior_list, delta_t,
#                            acf_list = list(rbf_acf, rbf_acf, rquad_acf),
#                            theta_prop_scale = theta_prop_scale,
#                            fix_gp_widthscale = fix_theta1,
#                            fix_param = fix_param)

if (do_sampling){
  #--------------- Start MCMC ---------------------
  cl <- makeCluster(n_chain)
  registerDoParallel(cl)
  cat("---------- Sampling started --------------\n")
  if (n_sam%%n_sam_batch_size != 0){
    stop("n_sam_batch_size must divide n_sam.")
  }
  sample_time <- system.time({
    all_samples <- foreach(i=1:n_chain, .packages = c("SuperGauss")) %dopar%{
      n_batch <- n_sam/n_sam_batch_size
      param_name_list <- c("X_sam","sig2_sam","beta_sam","alpha_sam","theta_sam")
      batch_id <- 1
      sam_batch_prev <- gpfa_gibbs_sampler(n_sam_batch_size, y_sim, init_param, 
                                           prior_list, delta_t, 
                                           acf_list = list(rbf_acf, rbf_acf, rquad_acf),
                                           theta_prop_scale = theta_prop_scale, 
                                           fix_gp_widthscale = fix_theta1,
                                           fix_param = fix_param)
      for (param_name in param_name_list){
        saveRDS(sam_batch_prev[[param_name]], batch_obj_name(param_name=param_name, chain_id=i, batch_id=batch_id))
      }
      for (batch_id in 2:n_batch){
        next_init_params <- list(X=sam_batch_prev$X_sam[n_sam_batch_size,,],
                                 beta=sam_batch_prev$beta_sam[n_sam_batch_size,,],
                                 alpha=sam_batch_prev$alpha_sam[n_sam_batch_size,],
                                 sig2=sam_batch_prev$sig2_sam[n_sam_batch_size,],
                                 theta=sam_batch_prev$theta_sam[n_sam_batch_size,])
        next_init_params$theta <- list(next_init_params$theta[1:2],
                                       next_init_params$theta[3:4],
                                       next_init_params$theta[5:7])
        sam_batch_prev <- gpfa_gibbs_sampler(n_sam_batch_size, y_sim, 
                                             next_init_params, 
                                             prior_list, delta_t, 
                                             acf_list = list(rbf_acf, rbf_acf, rquad_acf),
                                             theta_prop_scale = theta_prop_scale, 
                                             fix_gp_widthscale = fix_theta1,
                                             fix_param = fix_param)
        for (param_name in param_name_list){
          saveRDS(sam_batch_prev[[param_name]], batch_obj_name(param_name=param_name, chain_id=i, batch_id=batch_id))
        }
      }
    }
  })
  
  stopCluster(cl)
  cat("---------- Sampling completed --------------\n")
  print(sample_time)
}

if (do_sample_processing){
  source("gpfa_samples_processing.R")
}
