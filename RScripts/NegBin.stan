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

  vector[3] gamma;

  real<lower=0> phi1;
  real<lower=0> phi2;

  real<lower=0, upper=1> mix_weight;
}

transformed parameters {
  vector[N] mu1;
  vector[N] mu2;
  vector[N] phi1_vec;
  vector[N] phi2_vec;

  for (i in 1:N) {
    real season_sin = sin(2 * pi() * time[i] / 12);
    real season_cos = cos(2 * pi() * time[i] / 12);

    real linpred = alpha[region[i]] +
                   beta_time[region[i]] * time[i] +
                   beta_sin[region[i]] * season_sin +
                   beta_cos[region[i]] * season_cos +
                   theta_depth * depth[i] +
                   theta_lat   * lat[i] +
                   theta_lon   * lon[i] +
                   theta_lagmag * lag_mag[i];

    mu1[i] = exp(linpred);
    mu2[i] = exp(linpred + 0.75);  // heavier tail by boosting the mean

    phi1_vec[i] = exp(gamma[1] * nst[i] + gamma[2] * rms[i] + gamma[3] * clo[i]) * phi1;
    phi2_vec[i] = exp(gamma[1] * nst[i] + gamma[2] * rms[i] + gamma[3] * clo[i]) * phi2;
  }
}

model {
  // Priors
  alpha ~ normal(0, 1);
  beta_time ~ normal(0, 1);
  beta_sin ~ normal(0, 1);
  beta_cos ~ normal(0, 1);

  theta_depth ~ normal(0, 1);
  theta_lat ~ normal(0, 1);
  theta_lon ~ normal(0, 1);
  theta_lagmag ~ normal(0, 1);

  gamma ~ normal(0, 1);
  phi1 ~ exponential(1);
  phi2 ~ exponential(1);
  mix_weight ~ beta(2, 2);  // mildly regularized

  // Likelihood: mixture of NB1 and NB2
  for (i in 1:N) {
    target += log_mix(mix_weight,
      neg_binomial_2_lpmf(count[i] | mu1[i], phi1_vec[i]),
      neg_binomial_2_lpmf(count[i] | mu2[i], phi2_vec[i])
    );
  }
}

generated quantities {
  int y_rep[N];

  for (i in 1:N) {
    real component = bernoulli_rng(mix_weight);
    if (component == 1)
      y_rep[i] = neg_binomial_2_rng(mu1[i], phi1_vec[i]);
    else
      y_rep[i] = neg_binomial_2_rng(mu2[i], phi2_vec[i]);
  }
}