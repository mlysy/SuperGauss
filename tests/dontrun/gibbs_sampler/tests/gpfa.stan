data {
  int<lower=1> n_time;
  int<lower=1> n_proc;
  int<lower=1> n_factor;
  matrix[n_time, n_proc] Y;
  real timestamps[n_time];
  // prior list
  real sig_alpha[n_proc]; // shape for IGamma prior on var_1, ..., var_D
  real sig_beta[n_proc];  // rate for IGamma prior on var_1, ..., var_D
  vector[n_factor] psi_d;   // mean for Normal prior on beta_d
  matrix[n_factor, n_factor] S_d;  // var for Normal prior on beta_d
  vector[n_proc] nu; // mean for Normal prior on alpha
  vector[n_proc] tau; // var for Normal prior on alpha
}

parameters {
  matrix[n_time, n_factor] X;
  matrix[2, n_factor] theta;
  matrix[n_factor, n_proc] beta;
  row_vector[n_proc] alpha;
  vector[n_proc] sig2;
}

model {
  for (f in 1:n_factor){
    matrix[n_time, n_time] K_f = gp_exp_quad_cov(timestamps, theta[1, f], theta[2, f]);
    matrix[n_time, n_time] L_f = cholesky_decompose(K_f);
    X[,f] ~ multi_normal_cholesky(rep_vector(0, n_time), L_f);
  }
  
  for (t in 1:n_time){
    Y[t,] ~ multi_normal(X[t,]*beta + alpha, diag_matrix(sig2));
  }
  
  for (p in 1:n_proc){
    sig2[p] ~ inv_gamma(sig_alpha[p], sig_beta[p]);
    beta[,p] ~ multi_normal(psi_d, S_d);
    alpha[p] ~ normal(nu[p], sqrt(tau[p]));
  }
  
}
