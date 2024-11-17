require(MASS)
require(SuperGauss)
require(doParallel)
require(abind)
require(tidyverse)
require(ggpubr)
sampler_scripts <- list.files("gibbs_sampler", pattern="*.R$", 
                              full.names=TRUE, ignore.case=TRUE)
sapply(sampler_scripts, source, .GlobalEnv)

# Set what to save
save_figs <- TRUE
save_samples_as_rds <- TRUE

# Set MCMC parameters
n_chain <- 4
n_sam <- 3e4
n_warmup <- 1e4

# Set simulation parameters
fix_theta1 <- TRUE
set.seed(1234)
n_proc <- 30
n_factor <- 3
n_time <- 10000
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
par(mfrow=c(1,1))
plot(acf1, type="l", xlab="", ylab="acf", main="Black acf1 Blue acf2 Red acf3")
lines(acf2, col="blue")
lines(acf3, col="red")
x1 <- rnormtz(n=1, acf1)
x2 <- rnormtz(n=1, acf2)
x3 <- rnormtz(n=1, acf3)
x_sim <- cbind(x1, x2, x3)
# Plot x
xk_titles <- factor(c("Short memory X1", "Short memory X2", "Long memory X3"),
                    levels=c("Short memory X1", "Short memory X2", "Long memory X3"))
x_plot_df <- data.frame(latent_x=c(x1, x2, x3), 
                        time=rep(1:n_time, times=n_factor),
                        type=rep(xk_titles, each=n_time))
ggplot(data=x_plot_df)+
  geom_line(aes(x=time, y=latent_x))+facet_grid(rows=vars(type))+
  labs(x="Time", y="Simulated latent path")
if (save_figs) ggsave("ggpfa-simulated-x.pdf", width=10, height=6)

# Simulate observed data
y_sim <- x_sim %*% beta + 
  matrix(rep(alpha, n_time), nrow=n_time, ncol=n_proc, byrow = TRUE) +
  mvrnorm(n=n_time, mu=rep(0, n_proc), Sigma=diag(Sigma_diag))
cat("Number of simulated Y_i smaller than 0 is:\n",
    sum(y_sim < 0), "\n")
# Define priors and random initial values
prior_list <- list(sig_alpha=rep(1, n_proc),  # length D
                   sig_beta=rep(1, n_proc),   # length D
                   psi_d=rep(1, n_factor),      # length K
                   S_d=diag(rep(25,n_factor)),
                   nu=rep(0,n_proc),
                   tau=rep(100,n_proc),
                   ell_alpha=rep(NA, 3),
                   ell_beta=rep(NA, 3),
                   theta_additional_prior=function(theta){
                     # Only apply this prior on the rf kernel d parameter
                     if (length(theta)==3){
                       dunif(theta[3], -0.5, 0.5, log=TRUE)
                     } else{
                       0
                     }
                   })

# # Good initial values: true values
# init_param <- list(X=x_sim, 
#                    beta=beta, 
#                    alpha=alpha,
#                    sig2=Sigma_diag, 
#                    theta=theta)

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

# Better initial values: assuming we can get some noisy point estimates close to the truth
init_param <- list(X=x_sim+matrix(rnorm(length(x_sim),sd=0.05),nrow(x_sim),ncol(x_sim)), 
                   beta=beta+matrix(rnorm(length(beta),sd=1),nrow(beta),ncol(beta)), 
                   alpha=alpha+rnorm(length(alpha),sd=1),
                   sig2=Sigma_diag+rnorm(length(Sigma_diag),sd=0.05), 
                   theta=lapply(theta, function(e){e+c(0, rnorm(length(e)-1, sd=0.01))}))
# Tune theta_prop_scale
theta_prop_scale <- list(c(0, 0.001),  # first element set to 0 because fixing widthscale
                         c(0, 0.001), 
                         c(0, 0.004, 0.01))
# test <- gpfa_gibbs_sampler(300, y_sim, init_param, prior_list, delta_t,
#                            acf_list = list(rbf_acf, rbf_acf, rquad_acf),
#                            theta_prop_scale = theta_prop_scale,
#                            fix_gp_widthscale = fix_theta1)

#--------------- Start MCMC ---------------------
cl <- makeCluster(n_chain)
registerDoParallel(cl)
cat("---------- Sampling started --------------\n")
sample_time <- system.time({
  all_samples <- foreach(i=1:n_chain, .packages = c("SuperGauss")) %dopar%{
    gpfa_gibbs_sampler(n_sam, y_sim, init_param, prior_list, delta_t, 
                       acf_list = list(rbf_acf, rbf_acf, rquad_acf),
                       theta_prop_scale = theta_prop_scale, fix_gp_widthscale = fix_theta1)
  }
})
stopCluster(cl)
cat("---------- Sampling completed --------------\n")
print(sample_time)

#------------ Process the samples ----------------------
beta_sam <- abind(lapply(all_samples, function(chain) chain$beta_sam[-(1:n_warmup),,]), along=1)
theta_sam <- abind(lapply(all_samples, function(chain) chain$theta_sam[-(1:n_warmup),]), along=1)
sig2_sam <- abind(lapply(all_samples, function(chain) chain$sig2_sam[-(1:n_warmup),]), along=1)
alpha_sam <- abind(lapply(all_samples, function(chain) chain$alpha_sam[-(1:n_warmup),]), along=1)
X_sam <- abind(lapply(all_samples, function(chain) chain$X_sam[-(1:n_warmup),,]), along=1)
X_mean <- apply(X_sam, c(2,3), mean)
X_ci <- apply(X_sam, c(2,3), quantile, probs=c(0.025, 0.975))
if (fix_theta1){
  theta_sam <- theta_sam[, -c(1,3,5)]
}
sam_list <- list(beta_sam=beta_sam, theta_sam=theta_sam, sig2_sam=sig2_sam,
                 alpha_sam=alpha_sam, X_mean=X_mean, 
                 X_95ci_lo=X_ci["2.5%",,], X_95ci_hi=X_ci["97.5%",,])
if (save_samples_as_rds) {
  saveRDS(sam_list, "gaussian-gpfa-samples.rds")
  for (k in seq_len(n_factor)){
    saveRDS(X_sam[,,k], paste0("gaussian-gpfa-X-samples", k, ".rds"))
  }
  cat("------------- Samples saved ------------\n")
}
rm("beta_sam", "theta_sam", "sig2_sam", "alpha_sam", "X_sam")

#----------- Plots ----------------
#sam_list <- readRDS("gaussian-gpfa-samples.rds")
# Estimated posterior means
beta_mean <- apply(sam_list$beta_sam, c(2,3), mean)
beta_sd <- apply(sam_list$beta_sam, c(2,3), sd)
beta_ci_low <- apply(sam_list$beta_sam, c(2,3), quantile, probs=0.025)
beta_ci_hi <- apply(sam_list$beta_sam, c(2,3), quantile, probs=0.975)
sig2_mean <- apply(sam_list$sig2_sam, 2, mean)
sig2_sd <- apply(sam_list$sig2_sam, 2, sd)
sig2_ci_low <- apply(sam_list$sig2_sam, 2, quantile, probs=0.025)
sig2_ci_hi <- apply(sam_list$sig2_sam, 2, quantile, probs=0.975)
alpha_mean <- apply(sam_list$alpha_sam, 2, mean)
alpha_sd <- apply(sam_list$alpha_sam, 2, sd)
alpha_ci_low <- apply(sam_list$alpha_sam, 2, quantile, probs=0.025)
alpha_ci_hi <- apply(sam_list$alpha_sam, 2, quantile, probs=0.975)
theta_mean <- apply(sam_list$theta_sam, 2, mean)
theta_sd <- apply(sam_list$theta_sam, 2, sd)
# Estimated posterior means
Y_pred_mean <- sam_list$X_mean %*% beta_mean + 
  matrix(rep(alpha_mean, n_time), nrow=n_time, ncol=n_proc, byrow=TRUE)

data.frame(rbind(sig2_mean, sig2_sd, Sigma_diag), 
           row.names = c("Est_sig2", "Est_sig2_SD", "True_sig2"))
data.frame(rbind(alpha_mean, alpha_sd, alpha), 
           row.names = c("Est_alpha", "Est_alpha_SD", "True_alpha"))
data.frame(rbind(beta_mean, beta_sd, beta), 
           row.names = c(paste0("Est_beta", 1:n_factor), 
                         paste0("Est_beta_SD", 1:n_factor),
                         paste0("True_beta", 1:n_factor)))
true_theta <- theta
if (fix_theta1){
  true_theta <- lapply(true_theta, function(e){e[-1]})
}
data.frame(rbind(theta_mean, theta_sd, unlist(true_theta)), 
           row.names = c(paste0("Est_theta"), 
                         paste0("Est_theta_SD"),
                         paste0("True_theta")))

# Posterior density and trace plots for GP hyperparameters
n_sam_rm <- n_sam-n_warmup
iter_inds <- seq_len(n_sam_rm)
chain_inds <- split(iter_inds, ceiling(iter_inds/(n_sam_rm)))
theta_pos_plots <- list()
theta_trace_plots <- list()
for (k in 1:(n_factor)){
  theta_pos_plots[[k]] <- 
    ggplot(data=data.frame(theta=sam_list$theta_sam[,k]),
           aes(x=theta))+
    geom_density()+
    geom_vline(xintercept=theta[[k]][2], color="red", linetype="dashed")+
    labs(x="", y="", title=bquote(lambda[.(k)]))
  theta_trace_plots[[k]] <-
    ggplot(data=data.frame(iter=rep(seq_len(n_sam_rm), times=n_chain),
                           sam=sam_list$theta_sam[,k],
                           chain=unlist(lapply(seq_len(n_chain), 
                                               function(i)paste(rep("Chain", n_sam_rm), i)))))+
    geom_line(aes(x=iter, y=sam, col=chain))+
    labs(x="Iter", y="", title=bquote(lambda[.(k)]))
}
theta_pos_plots[[n_factor+1]] <- 
  ggplot(data=data.frame(theta=sam_list$theta_sam[,n_factor+1]),
         aes(x=theta))+
  geom_density()+
  geom_vline(xintercept=theta[[n_factor]][3], color="red", linetype="dashed")+
  labs(x="", y="", title=expression(d))
theta_trace_plots[[n_factor+1]] <-
  ggplot(data=data.frame(iter=rep(seq_len(n_sam_rm), times=n_chain),
                         sam=sam_list$theta_sam[,n_factor+1],
                         chain=unlist(lapply(seq_len(n_chain), 
                                             function(i)paste(rep("Chain", n_sam_rm), i)))))+
  geom_line(aes(x=iter, y=sam, col=chain))+
  labs(x="Iter", y="", title=expression(d))
ggarrange(plotlist=theta_pos_plots, ncol=4)
if (save_figs) ggsave("ggpfa-theta-pos-dist.pdf", width=10, height=2)
ggarrange(plotlist=theta_trace_plots, nrow=4)
if (save_figs) ggsave("ggpfa-theta-trace-plots.pdf", width=10, height=8)



alpha_df <- data.frame(est=alpha_mean, true=alpha,
                       ci_low=alpha_ci_low, ci_hi=alpha_ci_hi)
sig2_df <- data.frame(est=sig2_mean, true=Sigma_diag,
                      ci_low=sig2_ci_low, ci_hi=sig2_ci_hi)
beta_df <- data.frame(est=as.numeric(beta_mean), true=as.numeric(beta),
                      ci_low=as.numeric(beta_ci_low), ci_hi=as.numeric(beta_ci_hi))
ggplot(data=alpha_df, aes(x=true, y=est))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi), width=.1)+
  geom_abline(slope=1, linetype="dashed", col="red")+
  labs(x="True", y="Posterior mean", title=expression(alpha))
if (save_figs) ggsave("ggpfa-alpha-scatterplot.pdf", width=3, height=3)
ggplot(data=sig2_df, aes(x=true, y=est))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi), width=.1)+
  geom_abline(slope=1, linetype="dashed", col="red")+
  labs(x="True", y="Posterior mean", title=expression(sigma^2))
if (save_figs) ggsave("ggpfa-sig2-scatterplot.pdf", width=3, height=3)
ggplot(data=beta_df, aes(x=true, y=est))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_low, ymax=ci_hi), width=.01)+
  geom_abline(slope=1, linetype="dashed", col="red")+
  labs(x="True", y="Posterior mean", title=expression(beta))
if (save_figs) ggsave("ggpfa-beta-scatterplot.pdf", width=3, height=3)

n_time_to_show <- 500
x_est_plots <- list()
for (k in 1:n_factor){
  x_df <- data.frame(time=1:n_time_to_show, est=sam_list$X_mean[1:n_time_to_show,k], true=x_sim[1:n_time_to_show,k], 
                     ci_lo=sam_list$X_95ci_lo[1:n_time_to_show,k], ci_hi=sam_list$X_95ci_hi[1:n_time_to_show,k])
  x_est_plots[[k]] <- 
    ggplot(data=x_df)+
    geom_line(aes(x=time, y=true, color="True X"))+
    geom_line(aes(x=time, y=est, col="Pos. mean of X"))+
    geom_ribbon(aes(x=time, ymin=ci_lo, ymax=ci_hi), alpha=0.8, fill = "grey")+
    labs(x="Time", y="x", title=bquote(x[.(k)]))
}
ggarrange(plotlist=x_est_plots, nrow=3, common.legend = TRUE)
if (save_figs) ggsave("ggpfa-x-path-est.pdf", width=7, height=5)


# QQ plot for X
qq_plot <- function(zscore, title="", 
                    x_lab="Theoretical Quantile", 
                    y_lab="Sample Quantile") {
  tibble(z = zscore) %>%
    ggplot(aes(sample = z)) +
    stat_qq_line(linewidth = 1, color = "blue") +
    stat_qq() +
    ylab(y_lab) + xlab(x_lab)+
    ggtitle(title)
}
for (k in seq_len(n_factor)){
  x_sam_k <- readRDS(paste0("gaussian-gpfa-X-samples", k, ".rds"))
  qnorm_x <- rep(0, n_time)
  for (t in seq_len(n_time)){
    qnorm_x[t] <- qnorm(mean(x_sam_k[,t] <= x_sim[t, k]))
  }
  qq_plot(qnorm_x)
  if (save_figs) ggsave(paste0("ggpfa-x-qqplot-", k, ".pdf"), width=4, height=4)
}
