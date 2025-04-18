data {
  int<lower=1> N;                  // number of observations
  int<lower=1> R;                  // number of regions
  int<lower=1, upper=R> region[N]; // region index
  int<lower=0> y[N];               // count of earthquakes

  vector[N] time;       // TimeIndex (months since start)
  vector[N] depth;
  vector[N] lat;
  vector[N] lon;
  vector[N] lag_mag;

  vector[N] nst;
  vector[N] rms;
  vector[N] clo;
}


parameters {
  // Region-level coefficients
  vector[R] alpha;
  vector[R] beta_time;
  vector[R] beta_sin;
  vector[R] beta_cos;

  // Global coefficients for covariates
  real theta_depth;
  real theta_lat;
  real theta_lon;
  real theta_lagmag;

  // Overdispersion (phi) regression
  vector[3] gamma;

  // Hyperpriors for region effects
  real mu_alpha;
  real<lower=0> sigma_alpha;

  real mu_beta_time;
  real<lower=0> sigma_beta_time;

  real mu_beta_sin;
  real<lower=0> sigma_beta_sin;

  real mu_beta_cos;
  real<lower=0> sigma_beta_cos;
}

transformed parameters {
  vector[N] mu;
  vector[N] phi;

  for (i in 1:N) {
    real season_sin = sin(2 * pi() * time[i] / 12);
    real season_cos = cos(2 * pi() * time[i] / 12);

    mu[i] = exp(fmin(
      alpha[region[i]] +
      beta_time[region[i]] * time[i] +
      beta_sin[region[i]] * season_sin +
      beta_cos[region[i]] * season_cos +
      theta_depth * depth[i] +
      theta_lat * lat[i] +
      theta_lon * lon[i] +
      theta_lagmag * lag_mag[i],
      20  // prevents exp() overflow
      )
    );
    
    #if (is_nan(mu[i]) || mu[i] <= 0 || mu[i] > positive_infinity())
      #mu[i] = 1e-3;

    phi[i] = exp(
      fmin(gamma[1] * nst[i] +
      gamma[2] * rms[i] +
      gamma[3] * clo[i], 20)
    );
  }
}

model {
  // Hyperpriors
  mu_alpha ~ normal(0, 10);
  sigma_alpha ~ cauchy(0, 2);

  mu_beta_time ~ normal(0, 10);
  sigma_beta_time ~ cauchy(0, 2);

  mu_beta_sin ~ normal(0, 10);
  sigma_beta_sin ~ cauchy(0, 2);

  mu_beta_cos ~ normal(0, 10);
  sigma_beta_cos ~ cauchy(0, 2);

  // Region-level effects
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta_time ~ normal(mu_beta_time, sigma_beta_time);
  beta_sin ~ normal(mu_beta_sin, sigma_beta_sin);
  beta_cos ~ normal(mu_beta_cos, sigma_beta_cos);

  // Global covariates
  theta_depth ~ normal(0, 2);
  theta_lat ~ normal(0, 2);
  theta_lon ~ normal(0, 2);
  theta_lagmag ~ normal(0, 2);

  // Dispersion coefficients
  gamma ~ normal(0, 2);

  // Likelihood
  y ~ neg_binomial_2(mu, phi);
}

generated quantities {
  int y_rep[N];
  for (i in 1:N) {
    if (!is_nan(mu[i]) && !is_nan(phi[i]) &&
        mu[i] > 0 && mu[i] < 1e6 &&
        phi[i] > 0 && phi[i] < 1e6) {
      y_rep[i] = neg_binomial_2_rng(mu[i], phi[i]);
    } else {
      y_rep[i] = -1;
    }
  }
}

