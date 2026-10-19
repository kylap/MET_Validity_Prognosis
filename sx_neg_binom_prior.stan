//neg binomial model - symptoms
data {
  int<lower=0> N; 
  int<lower = 1> N_ids; 
  int<lower = 1, upper = 2> group[N]; 
  int<lower = 1, upper = 2> sex[N]; 
  int<lower = 1, upper = 7> timepoint[N]; 
  int<lower=1, upper=13> task[N];
  int<lower = 1, upper = 182> gr_sex_tp_ta[N]; 
  int<lower = 1> id [N]; 
  int<lower = 0, upper=162> sx[N]; 
}

parameters {
  vector[2] b_group; 
  vector[2] b_sex; 
  vector[7] b_time_point; 
  vector[13] b_task;
  vector[182] b_gr_sex_tp_ta;
  vector[N_ids] b_id_raw; 
  real mu_id;
  real<lower=0> sig_id; 
  real<lower = 0> phi; 
}

transformed parameters {
  vector[N_ids] b_id; 
  b_id = mu_id + sig_id * b_id_raw;
}

model {
  // Priors
  b_group[1] ~ normal(log(1),0.2); 
  b_group[2] ~ normal(log(4),0.3); 
  b_sex[1] ~ normal(log(1.5),0.2); 
  b_sex[2] ~ normal(log(2),0.2); 
  b_time_point[1] ~ normal(log(1),0.1); 
  b_time_point[2] ~ normal(log(4),0.3); 
  b_time_point[3] ~ normal(log(3),0.3); 
  b_time_point[4] ~ normal(log(2.5),0.2); 
  b_time_point[5] ~ normal(log(1),0.2); 
  b_time_point[6] ~ normal(log(1),0.2); 
  b_time_point[7] ~ normal(log(1),0.2); 
  b_task ~ normal(log(1.5),0.3); 
  
  b_gr_sex_tp_ta ~ normal(0,0.2); 
  
  mu_id ~ normal(0,0.3); 
  sig_id ~ normal(0.5,0.05); 
  b_id_raw ~ normal(0,1);
  
  phi ~ gamma(2,0.1);

  // Likelihood
  //vector[N] mu;
  //for (i in 1:N) {
  //  mu[i] = exp(b_group[group[i]] +
  //  b_sex[sex[i]] + 
  //  b_time_point[timepoint[i]] +
  //  b_task[task[i]] + 
  //  b_gr_sex_tp_ta[gr_sex_tp_ta[i]] +
  //  b_id[id[i]]
  //  );
  //}
  //sx ~ neg_binomial_2(mu, phi);
}

generated quantities {
  //vector[N] log_lik;
  //for (i in 1:N) {
  //  log_lik[i] = neg_binomial_2_lpmf(sx[i] | exp(b_group[group[i]] +
  //          b_sex[sex[i]] + 
  //          b_time_point[timepoint[i]] +
  //          b_task[task[i]] + 
  //          b_gr_sex_tp_ta[gr_sex_tp_ta[i]] +
  //          b_id[id[i]]
  //          ), phi);
  //}
}
