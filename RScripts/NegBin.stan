data {
  int<lower=1> N;
  int<lower=1> R;
  int<lower=1, upper=R> region[N];
  int<lower=0> count[N];

  vector[N] time;
  vector[N] depth;
  vector[N] lat;
  vector[N] lon;
  vector[N] lag_mag;
  vector[N] nst;
  vector[N] rms;
  vector[N] clo;
}

parameters {
  // Non-centered parameterization for hierarchical effects
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[R] alpha_raw;

  real mu_beta_time;
  real<lower=0> sigma_beta_time;
  vector[R] beta_time_raw;

  real mu_beta_sin;
  real<lower=0> sigma_beta_sin;
  vector[R] beta_sin_raw;

  real mu_beta_cos;
  real<lower=0> sigma_beta_cos;
  vector[R] beta_cos_raw;

  real theta_depth;
  real theta_lat;
  real theta_lon;
  real theta_lagmag;

  vector[3] gamma;

  real<lower=0> phi;
}

transformed parameters {
  vector[R] alpha = mu_alpha + sigma_alpha * alpha_raw;
  vector[R] beta_time = mu_beta_time + sigma_beta_time * beta_time_raw;
  vector[R] beta_sin = mu_beta_sin + sigma_beta_sin * beta_sin_raw;
  vector[R] beta_cos = mu_beta_cos + sigma_beta_cos * beta_cos_raw;

  vector[N] mu;
  vector[N] phi_vec;

  for (i in 1:N) {
    mu[i] = exp(
      alpha[region[i]] +
      beta_time[region[i]] * time[i] +
      beta_sin[region[i]] * sin(2 * pi() * time[i] / 12) +
      beta_cos[region[i]] * cos(2 * pi() * time[i] / 12) +
      theta_depth * depth[i] +
      theta_lat   * lat[i] +
      theta_lon   * lon[i] +
      theta_lagmag * lag_mag[i]
    );

    real log_phi = gamma[1] * nst[i] + gamma[2] * rms[i] + gamma[3] * clo[i];
    phi_vec[i] = exp(log_phi) * phi;
  }
}

model {
  // Heavy-tailed priors for flexibility
  alpha_raw ~ normal(0, 1);
  beta_time_raw ~ normal(0, 1);
  beta_sin_raw ~ normal(0, 1);
  beta_cos_raw ~ normal(0, 1);

  mu_alpha ~ student_t(3, 0, 2);
  sigma_alpha ~ student_t(3, 0, 2);
  mu_beta_time ~ student_t(3, 0, 2);
  sigma_beta_time ~ student_t(3, 0, 2);
  mu_beta_sin ~ student_t(3, 0, 2);
  sigma_beta_sin ~ student_t(3, 0, 2);
  mu_beta_cos ~ student_t(3, 0, 2);
  sigma_beta_cos ~ student_t(3, 0, 2);

  theta_depth ~ student_t(3, 0, 2);
  theta_lat ~ student_t(3, 0, 2);
  theta_lon ~ student_t(3, 0, 2);
  theta_lagmag ~ student_t(3, 0, 2);

  gamma ~ normal(0, 0.5);
  phi ~ exponential(1);

  count ~ neg_binomial_2(mu, phi_vec);
}

generated quantities {
  int y_rep[N];
  for (i in 1:N) {
    y_rep[i] = neg_binomial_2_rng(mu[i], phi_vec[i]);
  }
}
