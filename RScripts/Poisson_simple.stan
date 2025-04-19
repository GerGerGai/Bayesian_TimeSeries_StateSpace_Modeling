data {
  int<lower=1> N;                    // number of observations
  int<lower=0> count[N];             // response variable

  vector[N] time;        // continuous time index
  vector[N] depth;

  vector[N] lag_mag;
}

parameters {
  real alpha;                // intercept
  real<lower=-1,upper=1> phi_raw;
  real<lower=0> sigma;       // RW sd
  vector[N] z;               // non‑centred states
  real theta_depth;
  real theta_lagmag;

}

transformed parameters {
  real phi = phi_raw;        // already constrained
  vector[N] log_lambda;
  log_lambda[1] = alpha + z[1]*sigma;
  for (t in 2:N)
    log_lambda[t] = alpha + phi*(log_lambda[t-1]-alpha) + z[t]*sigma;
}

model {
  // Priors
  alpha        ~ normal(0, 1);
  phi_raw      ~ uniform(-1, 1);
  sigma        ~ cauchy(0, 1);     // half‑Cauchy automatically by <0>
  theta_depth   ~ normal(0, 0.5);
  theta_lagmag  ~ normal(0, 0.5);
  z            ~ normal(0, 1);


  // Likelihood
  for (t in 1:N)
    count[t] ~ poisson_log(log_lambda[t] +
                       theta_depth  * depth[t] +
                       theta_lagmag * lag_mag[t]);
}

generated quantities {
  int y_rep[N];

   for (t in 1:N)
    y_rep[t] = poisson_log_rng(log_lambda[t] +
                               theta_depth  * depth[t] +
                               theta_lagmag * lag_mag[t]);
}
