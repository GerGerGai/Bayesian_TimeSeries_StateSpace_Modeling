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
  vector[R] alpha;
  vector[R] beta_time;
  vector[R] beta_sin;
  vector[R] beta_cos;

  real theta_depth;
  real theta_lat;
  real theta_lon;
  real theta_lagmag;

  vector[3] gamma;     // for nst, rms, clo
  real<lower=0> phi_base;
}

transformed parameters {
  vector[N] mu;
  vector[N] phi;

  for (i in 1:N) {
    real season_sin = sin(2 * pi() * time[i] / 12);
    real season_cos = cos(2 * pi() * time[i] / 12);

    mu[i] = exp(
      alpha[region[i]] +
      beta_time[region[i]] * time[i] +
      beta_sin[region[i]] * season_sin +
      beta_cos[region[i]] * season_cos +
      theta_depth * depth[i] +
      theta_lat   * lat[i] +
      theta_lon   * lon[i] +
      theta_lagmag * lag_mag[i]
    );

    real log_phi = gamma[1] * nst[i] + gamma[2] * rms[i] + gamma[3] * clo[i];
    phi[i] = exp(log_phi) * phi_base;
  }
}

model {
  alpha ~ normal(0, 2);
  beta_time ~ normal(0, 2);
  beta_sin ~ normal(0, 2);
  beta_cos ~ normal(0, 2);

  theta_depth ~ normal(0, 2);
  theta_lat   ~ normal(0, 2);
  theta_lon   ~ normal(0, 2);
  theta_lagmag ~ normal(0, 2);

  gamma ~ normal(0, 0.5);      // regularized dispersion
  phi_base ~ exponential(1);  // base level of dispersion

  for (i in 1:N)
    count[i] ~ neg_binomial_2(mu[i], phi[i]);
}

generated quantities {
  vector[N] y_rep;
  for (i in 1:N) {
    y_rep[i] = neg_binomial_2_rng(mu[i], phi[i]);
  }
}

