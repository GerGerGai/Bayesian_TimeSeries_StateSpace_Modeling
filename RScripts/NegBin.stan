data {
  int<lower=1> N;                    // number of observations
  int<lower=1> R;                    // number of regions
  int<lower=1, upper=R> region[N];   // region ID for each observation
  int<lower=0> count[N];             // response variable

  vector[N] time;        // continuous time index
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

  vector[3] gamma;  // nst, rms, clo

  real<lower=0> phi_base;  // baseline overdispersion
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

    // Safety net
    if (is_nan(mu[i]) || mu[i] <= 0 || mu[i] > positive_infinity())
      mu[i] = 1e-3;

    phi[i] = exp(
      log(phi_base) +
      gamma[1] * nst[i] +
      gamma[2] * rms[i] +
      gamma[3] * clo[i]
    );
  }
}

model {
  // Priors
  alpha ~ normal(0, 2);
  beta_time ~ normal(0, 1);
  beta_sin ~ normal(0, 1);
  beta_cos ~ normal(0, 1);

  theta_depth ~ normal(0, 1);
  theta_lat ~ normal(0, 1);
  theta_lon ~ normal(0, 1);
  theta_lagmag ~ normal(0, 1);

  gamma ~ normal(0, 1);
  phi_base ~ exponential(1);

  // Likelihood
  count ~ neg_binomial_2(mu, phi);
}

generated quantities {
  int y_rep[N];

  for (i in 1:N) {
    if (is_nan(mu[i]) || mu[i] <= 0 || mu[i] > positive_infinity() ||
        is_nan(phi[i]) || phi[i] <= 0 || phi[i] > positive_infinity()) {
      y_rep[i] = -1;  // invalid prediction
    } else {
      y_rep[i] = neg_binomial_2_rng(mu[i], phi[i]);
    }
  }
}
