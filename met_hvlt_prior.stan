//MET HVLT - Prior
data {
  int<lower=0> N;    
  int<lower=1> N_ids;   
  int<lower=1, upper=2> group[N]; 
  int<lower=1, upper=2> sex[N]; 
  int<lower=1, upper=7> timepoint[N];
  int<lower=1, upper=14> gr_sex_tp[N]; 
  int<lower=1> id[N];   
  int<lower=0, upper=36> outcome[N]; 
  int<lower=1> M; 
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
}

transformed parameters {
  vector[14] b_gr_sex_tp;  
  vector[N_ids] b_id;  
  
  b_gr_sex_tp = mu_gr_sex_tp + b_gr_sex_tp_raw * sig_gr_sex_tp;
  b_id = mu_id + b_id_raw*sig_id;
}

model {
  // Priors
  b_group ~ normal(0.1,0.5); 
  b_sex[1] ~ normal(0.2198,0.5); 
  b_sex[2] ~ normal(0.3297,0.5);
  b_time_point[1] ~ normal(0.6594,0.5); 
  b_time_point[2] ~ normal(0.491,0.5);
  b_time_point[3] ~ normal(0.572,0.5); 
  b_time_point[4] ~ normal(0.6594,0.5); 
  b_time_point[5] ~ normal(0.6594,0.5); 
  b_time_point[6] ~ normal(0.6594,0.5); 
  b_time_point[7] ~ normal(0.6594,0.5); 
  mu_gr_sex_tp ~ normal(0, 0.5);
  sig_gr_sex_tp ~ exponential(1);
  b_gr_sex_tp_raw ~ normal(0, 1); 
  b_id_raw ~ normal(0, 1);    
  mu_id ~ normal(0, 0.5);     
  sig_id ~ exponential(1);   

  // Likelihood
  //for (i in 1:N) {
  //  real alpha = b_group[group[i]] +
  //             b_sex[sex[i]] +
  //             b_time_point[timepoint[i]] +
  //             b_gr_sex_tp[gr_sex_tp[i]] +
  //             b_id[id[i]];
  //  target += binomial_logit_lpmf(outcome[i] | M, alpha);
  //}
}

generated quantities {
  //vector[N] log_lik;
  //vector[N] pi;

  //for (i in 1:N) {
  //  real alpha = b_group[group[i]] +
  //             b_sex[sex[i]] +
  //             b_time_point[timepoint[i]] +
  //             b_gr_sex_tp[gr_sex_tp[i]] +
  //             b_id[id[i]];
  //  pi[i] = inv_logit(alpha); 
  //  log_lik[i] = binomial_logit_lpmf(outcome[i] | M, alpha);
  //}
}
