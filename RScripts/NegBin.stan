data {
  int<lower=0> N; // number of (month, region) groups
  int<lower=0> R; // number of unique regions
  int<lower=1, upper=R> region[N]; // region ID for each row (from 1 to R)
  int<lower=0> y[N]; // observed earthquake counts with magnitude >= 4
  vector[N] depth;
  vector[N] lat;
  vector[N] lon;
  vector[N] time;
  vector[N] lag_mag;
  vector[N] nst;
  vector[N] rms;
  vector[N] clo;
}


parameters {
  real beta_0, beta_1, beta_2, beta_3, beta_4, beta_5; // coefficients for predicting log mean
  real gamma_0, gamma_1, gamma_2, gamma_3; // coefficients for log dispersion
  vector[R] u; // random intercepts for each region
  real<lower=0> sigma_u;// sd of u
}

transformed parameters {
  vector[N] mu; // mean of NegBin
  vector[N] phi; // dispersion of NegBin
  for (i in 1:N) {
    mu[i] = exp(beta_0 + beta_1 * depth[i] + beta_2 * lat[i] + beta_3 * lon[i] +
                beta_4 * time[i] + beta_5 * lag_mag[i] + u[region[i]]);
    phi[i] = exp(gamma_0 + gamma_1 * nst[i] + gamma_2 * rms[i] + gamma_3 * clo[i]);
  }
}
model {
  
  beta_0 ~ normal(0, 10); 
  beta_1 ~ normal(0, 10);
  beta_2 ~ normal(0, 10); 
  beta_3 ~ normal(0, 10);
  beta_4 ~ normal(0, 10); 
  beta_5 ~ normal(0, 10);
  
  gamma_0 ~ normal(0, 10); 
  gamma_1 ~ normal(0, 10);
  gamma_2 ~ normal(0, 10); 
  gamma_3 ~ normal(0, 10);
  
  sigma_u ~ cauchy(0, 5); // cauchy distribution which has heavy tails and wide spreads.
  u ~ normal(0, sigma_u);
  
  for (i in 1:N) {
    y[i] ~ neg_binomial_2(mu[i], phi[i]);
  }
}
generated quantities {
  int y_rep[N]; // posterior predictive draws
  for (i in 1:N) {
    real safe_mu = fmin(mu[i], 1e6);       // cap mu
    real safe_phi = fmin(phi[i], 1e6);     // cap phi
    
    // Catch and fix bad values
    if (is_nan(safe_mu) || safe_mu > 1e6 || safe_mu < 0)
      safe_mu = 1e6;
    if (is_nan(safe_phi) || safe_phi > 1e6 || safe_phi < 0)
      safe_phi = 1e6;
      
    y_rep[i] = neg_binomial_2_rng(safe_mu, safe_phi);
  }
}