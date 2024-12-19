init_params_factanal <- function(y_sim, n_factor, acf_list, fixed_gp_var=0.05){
  n_proc <- ncol(y_sim)
  n_time <- nrow(y_sim)
  fa_results <- factanal(y_sim, n_factor, scores="regression")
  init_x <- fa_results$scores
  init_beta <- t(fa_results$loadings[1:n_proc,])
  init_sig2 <- fa_results$uniquenesses
  # Rearrange beta row order 
  beta_order <- c(which.max(init_beta[,1]), 
                  which.max(init_beta[,11]), 
                  which.max(init_beta[,21]))
  init_beta <- init_beta[beta_order, ]
  init_x <- init_x[, beta_order]
  
  # Scale variance of x based on the fixed GP variance
  for (d in 1:n_factor){
    scale_d <- sqrt(fixed_gp_var)/sd(init_x[,d])
    init_x[,d] <- init_x[,d]*scale_d
    init_beta[d,] <- init_beta[d,]/scale_d
  }
  
  # One draw from Gibbs sampler to initialize other parameters
  init_alpha <- rep(0, n_proc)
  for (d in 1:n_proc){
    init_alpha[d] <- cond_alpha_d(y_sim[,d], init_x, init_sig2[d], init_beta[,d], 
                                  prior_list$nu[d], prior_list$tau[d])$sample
  }
  XtX <- crossprod(init_x)
  for (d in 1:n_proc){
    init_beta[,d] <- cond_beta_d(y_sim[,d], init_x, XtX, init_sig2[d], init_alpha[d],
                                 prior_list$psi_d, prior_list$S_d)$sample
  }
  
  # Optimization to get GP hyperparameter init values
  obj1 <- function(lambda){
    -dnormtz(X=init_x[,1], 
             acf=acf_list[[1]](lag=(0:(n_time-1))*delta_t, 
                               var=fixed_gp_var, ell=lambda), 
             log=T, method="gschur")
  }
  obj2 <- function(lambda){
    -dnormtz(X=init_x[,2], 
             acf=acf_list[[2]](lag=(0:(n_time-1))*delta_t, 
                               var=fixed_gp_var, ell=lambda), 
             log=T, method="gschur")
  }
  obj3 <- function(lambda_and_d){
    -dnormtz(X=init_x[,3], 
             acf=acf_list[[3]]((0:(n_time-1))*delta_t, var=fixed_gp_var, 
                               k=lambda_and_d[1], d=lambda_and_d[2]), 
             log=T, method="gschur")
  }
  
  opt1 <- optimize(obj1, c(0, 0.2), maximum = FALSE)
  opt2 <- optimize(obj2, c(0, 0.2), maximum = FALSE)
  opt3 <- optim(c(0.1, 0.4), obj3)
  list(X=init_x,
       beta=init_beta,
       alpha=init_alpha,
       sig2=init_sig2,
       theta=list(c(fixed_gp_var, opt1$minimum),
                  c(fixed_gp_var, opt2$minimum),
                  c(fixed_gp_var, opt3$par)))
}


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

batch_obj_name <- function(chain_id, batch_id, param_name){
  paste0("gpfa_", param_name, "_sam_chain", chain_id, "_batch", batch_id, ".rds")
}

