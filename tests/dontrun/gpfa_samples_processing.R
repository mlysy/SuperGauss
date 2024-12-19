batch_ids <- 2:5 # the first batch is discarded as warm-up

#------------------ Process X samples -------------------
dir.create("X_samples_temp", showWarnings = FALSE)
for (chain_id in 1:n_chain){
  for (batch_id in batch_ids){
    X_sam_batch <- readRDS(batch_obj_name(chain_id=chain_id, batch_id=batch_id, param_name="X_sam"))
    for (k in seq_len(n_factor)){
      saveRDS(X_sam_batch[,,k], 
              file.path("X_samples_temp", 
                        batch_obj_name(chain_id=chain_id, batch_id=batch_id, param_name=paste0("X_sam", k))))
    }
  }
}

# QQ plot for X
for (k in seq_len(n_factor)){
  x_sam_k_list <- vector("list", n_chain*length(batch_ids))
  count <- 1
  for (chain_id in 1:n_chain){
    for (batch_id in batch_ids){
      x_sam_k_list[[count]] <- readRDS(file.path("X_samples_temp", 
                                          batch_obj_name(chain_id=chain_id, batch_id=batch_id, param_name=paste0("X_sam", k))))
      count <- count + 1
    }
  }
  x_sam_k <- do.call(rbind, x_sam_k_list)
  qnorm_x <- rep(0, n_time)
  for (t in seq_len(n_time)){
    qnorm_x[t] <- qnorm(mean(x_sam_k[,t] <= x_sim[t, k]))
  }
  qq_plot(qnorm_x)
  if (save_figs) ggsave(paste0("ggpfa-x-qqplot-", k, ".pdf"), width=4, height=4)
}
unlink("X_samples_temp", recursive = TRUE)

#----------------- Process theta samples -----------------
theta_sam_list <- vector("list", n_chain*length(batch_ids))
count <- 1
for (chain_id in 1:n_chain){
  for (batch_id in batch_ids){
    theta_batch <- readRDS(batch_obj_name(chain_id=chain_id, batch_id=batch_id, 
                                          param_name="theta_sam"))
    theta_sam_list[[count]] <- theta_batch[,c(2, 4, 6, 7)]
    count <- count+1
  }
}
theta_sam <- data.frame(do.call(rbind, theta_sam_list))
colnames(theta_sam) <- c("lam1", "lam2", "lam3", "d")
theta_sam$chain <- unlist(lapply(seq_len(n_chain), 
                                 function(i)paste(rep("Chain", n_sam_batch_size*length(batch_ids)), i)))
theta_sam$iter <- rep(seq_len(n_sam_batch_size*length(batch_ids)), times=n_chain)

# Posterior density and trace plots for GP hyperparameters
theta_pos_plots <- theta_trace_plots <- vector("list", n_factor+1)
for (k in 1:(n_factor)){
  theta_pos_plots[[k]] <- 
    ggplot(data=data.frame(theta=theta_sam[,k]),
           aes(x=theta))+
    geom_density()+
    geom_vline(xintercept=theta[[k]][2], color="red", linetype="dashed")+
    labs(x="", y="", title=bquote(lambda[.(k)]))
  theta_trace_plots[[k]] <-
    ggplot(data=data.frame(iter=theta_sam$iter,
                           sam=theta_sam[,k],
                           chain=theta_sam$chain))+
    geom_line(aes(x=iter, y=sam, col=chain))+
    labs(x="Iter", y="", title=bquote(lambda[.(k)]))
}
theta_pos_plots[[n_factor+1]] <-
  ggplot(data=data.frame(theta=theta_sam[,n_factor+1]),
         aes(x=theta))+
  geom_density()+
  geom_vline(xintercept=theta[[n_factor]][3], color="red", linetype="dashed")+
  labs(x="", y="", title=expression(d))
theta_trace_plots[[n_factor+1]] <-
  ggplot(data=data.frame(iter=theta_sam$iter,
                         sam=theta_sam[,n_factor+1],
                         chain=theta_sam$chain))+
  geom_line(aes(x=iter, y=sam, col=chain))+
  labs(x="Iter", y="", title=expression(d))
ggarrange(plotlist=theta_pos_plots, ncol=4)
if (save_figs) ggsave("ggpfa-theta-pos-dist.pdf", width=10, height=2)
ggarrange(plotlist=theta_trace_plots, nrow=4)
if (save_figs) ggsave("ggpfa-theta-trace-plots.pdf", width=10, height=8)


param_names <- c("sig2_sam", "beta_sam", "alpha_sam")
map_true_vals <- function(param){
  if (param=="sig2_sam"){
    return(Sigma_diag)
  } else if (param=="beta_sam"){
    return(beta)
  } else if (param=="alpha_sam"){
    return(alpha)
  }
}
for (param in param_names){
  batch_list <- vector("list", n_chain*length(batch_ids))
  count <- 1
  for (c_id in 1:n_chain){
    for (b_id in batch_ids){
      batch_list[[count]] <- readRDS(batch_obj_name(chain_id=chain_id, batch_id=batch_id, 
                                                    param_name=param))
      count <- count+1
    }
  }
  if (param == "beta_sam"){
    param_sam <- abind(batch_list, along=1)
    apply_dim <- c(2,3)
  } else{
    param_sam <- do.call(rbind, batch_list)
    apply_dim <- 2
  }
  computed_stats <- lapply(list(mean, sd, 
                                function(x) quantile(x, probs=0.025),
                                function(x) quantile(x, probs=0.975)),
                           function(stat) apply(param_sam, apply_dim, stat))
  names(computed_stats) <- c("mean", "sd", "ci_lo", "ci_hi") 
  plot_df <- data.frame(est=as.numeric(computed_stats[["mean"]]), 
                        ci_low=as.numeric(computed_stats[["ci_lo"]]), 
                        ci_hi=as.numeric(computed_stats[["ci_hi"]]),
                        true=as.numeric(map_true_vals(param)))
  ggplot(data=plot_df, aes(x=true, y=est))+
    geom_point()+
    geom_errorbar(aes(ymin=ci_low, ymax=ci_hi), width=.1)+
    geom_abline(slope=1, linetype="dashed", col="red")+
    labs(x="True", y="Posterior mean", title=param)
  if (save_figs) ggsave(paste0("ggpfa-", param, "-scatterplot.pdf"), width=3, height=3)
}
