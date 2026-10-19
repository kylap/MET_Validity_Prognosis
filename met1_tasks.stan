//MET 1 Tasks
data {
  int<lower=1> N;
  int<lower = 1> N_ids;
  int<lower = 1, upper = 2> group[N]; 
  int<lower = 1, upper = 2> sex[N];
  int<lower = 1, upper = 7> timepoint[N];
  int<lower=1, upper=3> task[N];
  int<lower = 1, upper = 42> gr_sex_tp_ta[N];
  int<lower = 1, upper = 6> task_order[N];
  int<lower = 1> id [N];
  vector[N] y;
}

parameters {
  real<lower=0> nu;
  vector[2] b_group;
  vector[2] b_sex;
  vector[7] b_time_point;
  vector[3] b_task;
  vector[42] b_gr_sex_tp_ta_raw;
  real mu_gr_sex_tp_ta;
  real<lower=0> sig_gr_sex_tp_ta;
  vector[6] b_task_order;
  vector[N_ids] b_id_raw;
  real mu_id;
  real<lower=0> sig_id;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N_ids] b_id;
  b_id = mu_id + sig_id * b_id_raw;
  
  vector[42] b_gr_sex_tp_ta;
  b_gr_sex_tp_ta = mu_gr_sex_tp_ta + b_gr_sex_tp_ta_raw * sig_gr_sex_tp_ta;
}

model {
  // Priors
  b_group[1] ~ normal(0.2,0.2); 
  b_group[2] ~ normal(-0.2,0.2);
  b_sex[1] ~ normal(-0.1,0.2);
  b_sex[2] ~ normal(0,0.2);
  b_time_point ~ normal(0,0.2);
  b_task[1] ~ normal(-0.3,0.2); 
  b_task[2] ~ normal(0.1,0.2); 
  b_task[3] ~ normal(0.2,0.2);
  
  sig_gr_sex_tp_ta ~ exponential(5);
  mu_gr_sex_tp_ta ~ normal(0,0.5);
  b_gr_sex_tp_ta_raw ~ normal(0,1);
  
  b_task_order ~ normal(0,0.3);
  b_id_raw ~ normal(0, 1);
  mu_id ~ normal(0, 0.5);
  sig_id ~ exponential(5); 
  sigma ~ exponential(1);
  nu ~ gamma(2, 0.1);
  
  // Likelihood
  for (i in 1:N) {
    y[i] ~ student_t(nu, 
                     b_group[group[i]] +
                     b_sex[sex[i]] + 
                     b_time_point[timepoint[i]] +
                     b_task[task[i]] + 
                     b_gr_sex_tp_ta[gr_sex_tp_ta[i]] +
                     b_task_order[task_order[i]] +
                     b_id[id[i]], 
                     sigma);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;

  for (i in 1:N) {mu[i] = b_group[group[i]] +
                     b_sex[sex[i]] + 
                     b_time_point[timepoint[i]] +
                     b_task[task[i]] + 
                     b_gr_sex_tp_ta[gr_sex_tp_ta[i]] +
                     b_task_order[task_order[i]] +
                     b_id[id[i]];
  }
  for (i in 1:N) {
    log_lik[i] = student_t_lpdf(y[i] | nu, mu[i], sigma);
  }
}
