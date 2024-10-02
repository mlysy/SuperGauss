require(MASS)
require(SuperGauss)
require(doParallel)
require(abind)
require(tidyverse)
require(ggpubr)
sampler_scripts <- list.files("gibbs_sampler", pattern="*.R$", 
                              full.names=TRUE, ignore.case=TRUE)
sapply(sampler_scripts, source, .GlobalEnv)
save_figs <- FALSE
fix_theta1 <- TRUE
set.seed(123)
n_proc <- 12
n_factor <- 3
n_time <- 2000
delta_t <- 0.1
Sigma_diag <- rnorm(n_proc, 0.5, 0.1)^2
proc_per_group <- n_proc/n_factor
beta <- rbind(c(runif(proc_per_group, 3, 5), runif(proc_per_group, 0, 1), runif(proc_per_group, 0, 1)),
              c(runif(proc_per_group, 0, 1), runif(proc_per_group, 3, 5), runif(proc_per_group, 0, 1)),
              c(runif(proc_per_group, 0, 1), runif(proc_per_group, 0, 1), runif(proc_per_group, 3, 5)))
alpha <- rnorm(n_proc, mean=7, sd=2)
theta <- list(c(0.05, 0.05),
              c(0.05, 0.1),
              c(0.05, 0.12, 0.1))

# Simulate factor processes
# x_sim <- sapply(1:n_factor,
#                 function(k){
#                   x_acf <- rbf_acf((0:(n_time-1))*delta_t,
#                                    var=theta[[k]][1], ell=theta[[k]][2])
#                   rnormtz(n = 1, x_acf)
#                 })
# Short memory latent process
acf1 <- rbf_acf((0:(n_time-1))*delta_t, var=theta[[1]][1], ell=theta[[1]][2])
acf2 <- rbf_acf((0:(n_time-1))*delta_t, var=theta[[2]][1], ell=theta[[2]][2])
#acf3 <- arfima_acf((0:(n_time-1))*delta_t, var=theta[[3]][1], d=theta[[3]][2])
acf3 <- rquad_acf((0:(n_time-1))*delta_t, var=theta[[3]][1], k=theta[[3]][2], alpha=theta[[3]][3])
par(mfrow=c(1,1))
plot(acf1, type="l", xlab="", ylab="acf", main="Black acf1 Blue acf2 Red acf3")
lines(acf2, col="blue")
lines(acf3, col="red")
x1 <- rnormtz(n=1, acf1)
x2 <- rnormtz(n=1, acf2)
x3 <- rnormtz(n=1, acf3)
x_sim <- cbind(x1, x2, x3)
par(mfrow=c(3,1))
xk_titles <- c("Short memory x1", "Short memory x2", "Long memory x3")
for (k in 1:n_factor) {
  plot(x_sim[,k],type="l", ylab=paste0("x", k), main=xk_titles[k])
}

# Simulate observed data
y_sim <- x_sim %*% beta + 
  matrix(rep(alpha, n_time), nrow=n_time, ncol=n_proc, byrow = TRUE) +
  mvrnorm(n=n_time, mu=rep(0, n_proc), Sigma=diag(Sigma_diag))
plot(y_sim[,1], type="l", ylab="y1")
plot(y_sim[,2], type="l", ylab="y2")
plot(y_sim[,7], type="l", ylab="y7")
range(y_sim)
# Define priors and random initial values
prior_list <- list(sig_alpha=rep(1, n_proc),  # length D
                   sig_beta=rep(1, n_proc),   # length D
                   psi_d=rep(1, n_factor),      # length K
                   S_d=diag(rep(25,n_factor)),
                   nu=rep(0,n_proc),
                   tau=rep(100,n_proc),
                   ell_alpha=rep(NA, 3),
                   ell_beta=rep(NA, 3))
init_param <- list(X=x_sim, 
                   beta=beta, 
                   alpha=alpha,
                   sig2=Sigma_diag, 
                   theta=theta)
# More random initial values for x, just to test how the sampler performs...
init_x <- matrix(rep(0, n_time*n_factor), nrow=n_time, ncol=n_factor)
for (i in 1:n_factor){
  init_x[,i] <- ksmooth(1:n_time, y_sim[,i], "normal", bandwidth = 20)$y-alpha[i]
}
init_param <- list(X=init_x,
                   beta=matrix(1, nrow=n_factor, ncol=n_proc),
                   alpha=rep(1, n_proc),
                   sig2=rep(1, n_proc),
                   theta=list(rep(0.05, 2), rep(0.05, 2), rep(0.05, 3)))
init_param <- list(X=x_sim+matrix(rnorm(length(x_sim),sd=0.05),nrow(x_sim),ncol(x_sim)), 
                   beta=beta+matrix(rnorm(length(beta),sd=1),nrow(beta),ncol(beta)), 
                   alpha=alpha+rnorm(length(alpha),sd=1),
                   sig2=Sigma_diag+rnorm(length(Sigma_diag),sd=0.05), 
                   theta=lapply(theta, function(e){e+c(0, rnorm(length(e)-1, sd=0.01))}))
# Tune theta_prop_scale to get an acceptance rate around 0.4
theta_prop_scale <- list(0.001, 0.001, 0.005)
test <- gpfa_gibbs_sampler(100, y_sim, init_param, prior_list, delta_t, 
                           acf_list = list(rbf_acf, rbf_acf, rquad_acf),
                           theta_prop_scale = theta_prop_scale, 
                           fix_gp_widthscale = fix_theta1)

# Start MCMC
n_chain <- 2
n_sam <- 500
n_warmup <- 10
cl <- makeCluster(n_chain)
registerDoParallel(cl)
sample_time <- system.time({
  all_samples <- foreach(i=1:n_chain, .packages = c("SuperGauss")) %dopar%{
    gpfa_gibbs_sampler(n_sam, y_sim, init_param, prior_list, delta_t, 
                       acf_list = list(rbf_acf, rbf_acf, rquad_acf),
                       theta_prop_scale = theta_prop_scale, fix_gp_widthscale = fix_theta1)
  }
})
stopCluster(cl)

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
rm("beta_sam", "theta_sam", "sig2_sam", "alpha_sam", "X_sam")
#saveRDS(sam_list, "gaussian-gpfa-samples.rds")

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

theta_pos_plots <- list()
theta_trace_plots <- list()
for (k in 1:(n_factor)){
  theta_pos_plots[[k]] <- 
    ggplot(data=data.frame(theta=sam_list$theta_sam[,k]),
           aes(x=theta))+
    geom_density()+
    geom_vline(xintercept=theta[[k]][2], color="red", linetype="dashed")+
    labs(x="", y="", title=bquote(theta[.(k)]^ell))
  tot_sam <- length(sam_list$theta_sam[,k])
  theta_trace_plots[[k]] <-
    ggplot(data=data.frame(iter=seq_len(tot_sam/2),
                           chain1=sam_list$theta_sam[1:(tot_sam/2),k],
                           chain2=sam_list$theta_sam[(tot_sam/2+1):tot_sam,k]))+
    geom_line(aes(x=iter, y=chain1, col="Chain1"))+
    geom_line(aes(x=iter, y=chain2, col="Chain2"))+
    labs(x="Iter", y="", title=bquote(theta[.(k)]^ell))
}
theta_pos_plots[[n_factor+1]] <- 
  ggplot(data=data.frame(theta=sam_list$theta_sam[,n_factor+1]),
         aes(x=theta))+
  geom_density()+
  geom_vline(xintercept=theta[[n_factor]][3], color="red", linetype="dashed")+
  labs(x="", y="", title=expression(alpha))
theta_trace_plots[[n_factor+1]] <-
  ggplot(data=data.frame(iter=seq_len(tot_sam/2),
                         chain1=sam_list$theta_sam[1:(tot_sam/2),n_factor+1],
                         chain2=sam_list$theta_sam[(tot_sam/2+1):tot_sam,n_factor+1]))+
  geom_line(aes(x=iter, y=chain1, col="Chain1"))+
  geom_line(aes(x=iter, y=chain2, col="Chain2"))+
  labs(x="Iter", y="", title=expression(alpha))
ggarrange(plotlist=theta_pos_plots, ncol=4)
if (save_figs) ggsave("ggpfa-theta-pos-dist.pdf", width=10, height=2)
ggarrange(plotlist=theta_trace_plots, nrow=4)
if (save_figs) ggsave("ggpfa-theta-trace-plots.pdf", width=2, height=10)



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

x_est_plots <- list()
for (k in 1:n_factor){
  x_df <- data.frame(time=1:n_time, est=sam_list$X_mean[,k], true=x_sim[,k], 
                     ci_lo=sam_list$X_95ci_lo[,k], ci_hi=sam_list$X_95ci_hi[,k])
  x_est_plots[[k]] <- 
    ggplot(data=x_df)+
    geom_line(aes(x=time, y=true, color="True X"), alpha=0.8)+
    geom_line(aes(x=time, y=est, col="Pos. mean of X"), alpha=0.8)+
    geom_ribbon(aes(x=time, ymin=ci_lo, ymax=ci_hi, linetype=NA), alpha=0.4)+
    labs(x="Time", y="x", title=bquote(x[.(k)]))
}
ggarrange(plotlist=x_est_plots, nrow=3, common.legend = TRUE)
if (save_figs) ggsave("ggpfa-x-path-est.pdf", width=7, height=5)


#--------------------- Test against Stan samples ----------------------------
# require(rstan)
# require(dplyr)
# n_sam <- 5000
# n_warmup <- 1000
# n_chain <- 4
# mod <- stan_model(file="gpfa.stan")
# options(mc.cores=n_chain)
# stan_sam <- sampling(mod, chains=n_chain, iter=n_sam, warmup=n_warmup,
#                      data=c(list(n_time=n_time, n_proc=n_proc, n_factor=n_factor, 
#                                  Y=y_sim, 
#                                  timestamps=(1:n_time)*delta_t), prior_list),
#                      init=function() init_param)
# saveRDS(stan_sam, "stan_ggpfa.rds")
# all_samples <- data.frame(as.matrix(stan_sam))
# fill_ind <- ((nrow(all_samples)+1):(nrow(all_samples)+n_sam-n_warmup))
# all_samples$Method <- "Stan"
# all_samples <- all_samples %>% 
#   select(-lp__) %>%
#   mutate(Method="Stan")
# all_samples[fill_ind, ] <- NA
# all_samples$Method[fill_ind] <- "Gibbs"
# all_samples <- all_samples %>% select(-starts_with("X."))
# # Add samples from the Gibbs sampler
# for (ff in 1:n_factor){
#   for (pp in 1:n_proc){
#     all_samples[fill_ind, paste("beta", ff, pp, "" ,sep=".")] <- sam_list$beta_sam[(n_warmup+1):n_sam, ff, pp]
#   }
#   all_samples[fill_ind, paste("theta", 1, ff, "", sep=".")] <- sam_list$theta_sam[(n_warmup+1):n_sam, 1, ff]
#   all_samples[fill_ind, paste("theta", 2, ff, "", sep=".")] <- sam_list$theta_sam[(n_warmup+1):n_sam, 2, ff]
# }
# for (pp in 1:n_proc){
#   all_samples[fill_ind, paste("sig2", pp, "", sep=".")] <- sam_list$sig2_sam[(n_warmup+1):n_sam, pp]
# }
# 
# # Compare posterior dist in plots
# vars_to_plot <- c("theta.2.1.", "theta.2.2.", "theta.2.3.", colnames(all_samples)[7:66])
# true_vals <- c(as.vector(theta[2,]), as.vector(beta), alpha, Sigma_diag)
# test_plot_list <- vector("list", length(vars_to_plot))
# for (ii in 1:length(vars_to_plot)){
#   test_plot_list[[ii]] <- ggplot(data=all_samples, 
#                                  aes_string(x=vars_to_plot[ii], col="Method"))+
#     geom_density()+geom_vline(xintercept = true_vals[ii],linetype="dotted",color="green")+
#     xlab(vars_to_plot[ii])
#   print(test_plot_list[[ii]])
# }
