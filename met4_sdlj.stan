//Step-Down + Choice Lateral Jump - Binomial Regression
data {
  int<lower=0> N;   
  int<lower=1> N_ids;  
  int<lower=1, upper=2> group[N];
  int<lower=1, upper=2> sex[N]; 
  int<lower=1, upper=7> timepoint[N]; 
  int<lower=1, upper=14> gr_sex_tp[N]; 
  int<lower=1> id[N]; 
  int<lower=0, upper=12> outcome[N];
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
  
  for (i in 1:N_ids) {
    b_id[i] = mu_id + sig_id * b_id_raw[i];
  }
}

model {
  // Priors
  b_group ~ normal(logit(0.8),1);      
  b_sex ~ normal(logit(0.8),1);       
  b_time_point ~ normal(logit(0.8),1); 
  mu_gr_sex_tp ~ normal(0, 0.5);
  sig_gr_sex_tp ~ exponential(1);
  b_gr_sex_tp_raw ~ normal(0, 1); 
  b_id_raw ~ normal(0, 1);     
  mu_id ~ normal(0, 0.5);      
  sig_id ~ exponential(1);     

  // Likelihood
  for (i in 1:N) {
    real alpha = b_group[group[i]] +
               b_sex[sex[i]] +
               b_time_point[timepoint[i]] +
               b_gr_sex_tp[gr_sex_tp[i]] +
               b_id[id[i]];
    target += binomial_logit_lpmf(outcome[i] | M, alpha);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] pi;

  for (i in 1:N) {
    real alpha = b_group[group[i]] +
               b_sex[sex[i]] +
               b_time_point[timepoint[i]] +
               b_gr_sex_tp[gr_sex_tp[i]] +
               b_id[id[i]];
    pi[i] = inv_logit(alpha);  // Probability of success
    log_lik[i] = binomial_logit_lpmf(outcome[i] | M, alpha);
  }
}
