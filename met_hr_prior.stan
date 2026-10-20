//Prior: Looking at concussion vs. healthy between all the time points - MET Ave HR
data {
  int<lower=0> N; 
  int<lower = 1> N_ids; 
  int<lower = 1, upper = 2> group[N]; 
  int<lower = 1, upper = 2> sex[N]; 
  int<lower = 1, upper = 7> timepoint[N]; 
  int<lower = 1, upper = 5> stage[N]; 
  int<lower = 1, upper = 70> gr_sex_tp_st[N]; 
  int<lower = 1> id [N];
  vector[N] hr; 
}

parameters {
  real<lower=0> nu; // degrees of freedom for Student-t distribution
  vector[2] b_group; //coefficients for group
  vector[2] b_sex; // coefficients for sex
  vector[7] b_time_point; // coefficients for time points
  vector[5] b_stage; //coefficients for MET stages
  real mu_gr_sex_tp_st;
  real<lower=0> sig_gr_sex_tp_st; // adaptive scale for interactions
  vector[70] b_gr_sex_tp_st_raw; //coefficients for group, sex, tp interaction
  vector[N_ids] b_id_raw; // raw participant random effects
  real mu_id; // mean of random effects for participants
  real<lower=0> sig_id; // standard deviation of participant random effects
  real<lower=0> sigma; // residual error
}

transformed parameters {
  vector[N_ids] b_id; //participant ids
  b_id = mu_id + sig_id * b_id_raw;
    
  vector[70] b_gr_sex_tp_st; // Interaction effects with non-centered parameterization
  b_gr_sex_tp_st = mu_gr_sex_tp_st + b_gr_sex_tp_st_raw * sig_gr_sex_tp_st;
}

model {
  // Priors
  b_group ~ normal(0,0.3); 
  b_sex[1] ~ normal(-0.1,0.3); 
  b_sex[2] ~ normal(0.1,0.3); 
  b_time_point ~ normal(0,0.3); 
  b_stage[1] ~ normal(-1.2,0.1); 
  b_stage[2] ~ normal(0.1,0.2); 
  b_stage[3] ~ normal(0.1,0.2); 
  b_stage[4] ~ normal(0.1,0.2); 
  b_stage[5] ~ normal(1.4,0.2); 
  
  sig_gr_sex_tp_st ~ exponential(1);
  mu_gr_sex_tp_st ~ normal(0,0.2);
  b_gr_sex_tp_st_raw ~ normal(0,1);
  b_id_raw ~ normal(0, 1); 
  mu_id ~ normal(0, 0.5); 
  sig_id ~ exponential(1); 
  sigma ~ exponential(1); 
  nu ~ gamma(2, 0.1); 
  
  // Likelihood
  //for (i in 1:N) {
  //  hr[i] ~ student_t(nu, 
  //                   b_group[group[i]] + 
  //                   b_time_point[timepoint[i]] + 
  //                   b_sex[sex[i]] + 
  //                   b_stage[stage[i]] +
  //                   b_gr_sex_tp_st[gr_sex_tp_st[i]] +
  //                   b_id[id[i]],  
  //                   sigma);
  //}
}

generated quantities {
  //vector[N] log_lik;
  //vector[N] mu;

  //for (i in 1:N) {
  //  mu[i] = b_group[group[i]] + 
  //          b_time_point[timepoint[i]] + 
  //          b_sex[sex[i]] + 
  //          b_stage[stage[i]] +
  //          b_gr_sex_tp_st[gr_sex_tp_st[i]] +
  //          b_id[id[i]];
  //}
  //for (i in 1:N) {
  //  log_lik[i] = student_t_lpdf(hr[i] | nu, mu[i], sigma);
  //}
}
