//neg binomial model for MET Stage 3 COWAT Words
data {
  int<lower=0> N;
  int<lower = 1> N_ids; 
  int<lower = 1, upper = 2> group[N]; 
  int<lower = 1, upper = 2> sex[N]; 
  int<lower = 1, upper = 7> timepoint[N]; 
  int<lower = 1, upper = 14> gr_sex_tp[N]; 
  int<lower = 1> id [N]; 
  int<lower = 0> outcome[N]; 
}

parameters {
  vector[2] b_group; 
  vector[2] b_sex; 
  vector[7] b_time_point; 
  vector[14] b_gr_sex_tp_raw; 
  real mu_gr_sex_tp;
  real<lower=0> sig_gr_sex_tp;
  vector[N_ids] b_id_raw; 
  real mu_id; 
  real<lower=0> sig_id; 
  real<lower = 0> phi; 
}

transformed parameters {
  vector[N_ids] b_id; 
  b_id = mu_id + sig_id * b_id_raw;
  
  vector[14] b_gr_sex_tp; 
  b_gr_sex_tp = mu_gr_sex_tp + b_gr_sex_tp_raw * sig_gr_sex_tp;
}

model {
  // Priors
  b_group ~ normal(log(2),0.1); 
  b_sex[1] ~ normal(log(2.3),0.1); 
  b_sex[2] ~ normal(log(2.5),0.1); 
  b_time_point ~ normal(log(2),0.1); 
  
  mu_gr_sex_tp ~ normal(log(1),0.1);
  sig_gr_sex_tp ~ exponential(5);
  b_gr_sex_tp_raw ~ normal(0,1); 
  
  b_id_raw ~ normal(0, 1); 
  mu_id ~ normal(0,0.1); 
  sig_id ~ exponential(5); 
  phi ~ gamma(10,1);

  // Likelihood
  for (i in 1:N) {
    target += neg_binomial_2_lpmf(outcome[i] | exp(
      b_group[group[i]] + 
      b_sex[sex[i]] + 
      b_time_point[timepoint[i]] + 
      b_gr_sex_tp[gr_sex_tp[i]] +
      b_id[id[i]]),
      phi);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;

  for (i in 1:N) {
    mu[i] =   
      b_group[group[i]] + 
      b_sex[sex[i]] + 
      b_time_point[timepoint[i]] + 
      b_gr_sex_tp[gr_sex_tp[i]] +
      b_id[id[i]];
  }
  
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_lpmf(outcome[i] | exp(mu[i]),phi);
  }
}
