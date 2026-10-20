// NBM Survival model

data {
  int<lower=0> N;
  int<lower=0> y[N];             // observed recovery days (or censoring time)
  int<lower=0, upper=1> censored[N]; // 1 = censored, 0 = recovered
  int<lower=1, upper=2> group[N];
}

parameters {
vector[2] alpha;                    // intercept on log scale
  vector<lower=0>[2] phi;             // overdispersion parameter
}

transformed parameters {
  vector<lower=0>[N] mu;
  for (n in 1:N) {
    mu[n] = exp(alpha[group[n]]);
  }
}


model {
  // Prior only for now
  alpha[1] ~ normal(log(35), 0.5);    //centered on 35 days - fail
  alpha[2] ~ normal(log(14), 0.5);    //centered on 21 days - pass
  phi ~ gamma(2, 2);             // prior with moderate dispersion

  // Likelihood
  for (n in 1:N) {
    int g = group[n];
    if (censored[n] == 0)
      target += neg_binomial_2_lpmf(y[n] | mu[n], phi[g]);
    else
      target += neg_binomial_2_lccdf(y[n] | mu[n], phi[g]);
  }
}

generated quantities {
  vector[201] R_group1;
  vector[201] R_group2;
  vector[N] log_lik;
  
  // Simulated recovery times per person
  int y_sim[N];

  for (t in 0:200) {
    R_group1[t + 1] = 1 - neg_binomial_2_cdf(t, exp(alpha[1]), phi[1]);
    R_group2[t + 1] = 1 - neg_binomial_2_cdf(t, exp(alpha[2]), phi[2]);
  }
  
  for (n in 1:N) {
    y_sim[n] = neg_binomial_2_rng(exp(alpha[group[n]]), phi[group[n]]);
  }
  
  for (n in 1:N) {
    int g = group[n];
    if (censored[n] == 0) {
      log_lik[n] = neg_binomial_2_lpmf(y[n] | mu[n], phi[g]);
    } else {
      log_lik[n] = neg_binomial_2_lccdf(y[n] | mu[n], phi[g]);
    }
  }
}


