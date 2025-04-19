library(dplyr)
library(lubridate)
library(rstan)
library(tidyr)
library(bayesplot)
library(scoringRules)
library(posterior)   
#install.packages("V8")


## Data Preparation
df <- read.csv("cleaned_earthquake_data.csv")
df$datetime <- as.POSIXct(paste(df$Date, df$Time), format="%Y/%m/%d %H:%M:%OS", tz="UTC") # combine dates and times


df <- df %>%
  filter(!is.na(datetime) & Mag >= 3.5) 

df <- df %>%
  mutate(LagMag = lag(Mag)) %>%
  ungroup()



agg_df <- df %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarise(
    Count = n(),
    Depth = mean(Depth, na.rm = TRUE),
    LagMag = mean(LagMag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  drop_na()

agg_df <- agg_df %>%
  mutate(across(c(Depth, LagMag), scale)) %>%
  mutate(across(everything(), ~ifelse(is.finite(.), ., 1e-3)))





###########################################################################
####### SECTION: model 1: Poisson model ######################

H = 12 # forecast horizon
N_total = nrow(agg_df)
N_train = N_total - H


stan_data_1 <- list(
  N = N_train,
  count = agg_df$Count[1:N_train],
  depth = as.vector(agg_df$Depth)[1:N_train],
  lag_mag = as.vector(agg_df$LagMag)[1:N_train],
  
  H = H,
  depth_future = as.vector(agg_df$Depth)[(N_train+1):N_total],
  lag_mag_future = as.vector(agg_df$Depth)[(N_train+1):N_total]
)


init_fn_1 <- function() {
  R <- stan_data_1$R
  
  list(
    alpha = 0,
    phi_raw = 0,
    sigma = 0,
    mu_z = 0,
    sigma_z = 1,          
    theta_depth = 0,
    theta_lagmag = 0
    
  )
}

fit1 <- stan(
  file = "RScripts/Poisson_simple.stan",
  data = stan_data_1,
  iter = 2000,
  chains = 4,
  seed = 42,
  control = list(adapt_delta = 0.8, max_treedepth = 10),
  init = "random"
)

########  model1 posterior check:

# 0. fast slow mixing:
mcmc_trace(fit1, pars = "alpha") 
mcmc_trace(fit1, pars = "phi_raw") 
mcmc_trace(fit1, pars = "sigma")


# 1. convergence: Rhat and neff
# aim: Rhat<1.01, >400
summ = summary(fit1)$summary
bad_rhat = sum(summ[,"Rhat"] > 1.01)
bad_neff  =  sum(summ[,"n_eff"] < 400)

ratios1 <- neff_ratio(fit1)
print(ratios1)

cat("No. of bad Rhat:", bad_rhat, "\n")
cat("No. of bad n_eff:", bad_neff, "\n")

# 2. sampler dignostic: divergent transitions and tree depth
# aim: 0 div, 0 hit max depth

sampler_pars = rstan::get_sampler_params(fit1, inc_warmup = FALSE)
sampler_mat  = do.call(rbind, sampler_pars)

div_total = sum(sampler_mat[, "divergent__"])
cat("Total divergences:", div_total, "\n")

max_td = max(sampler_mat[, "treedepth__"])
cat("Maximum treedepth reached:", max_td, "\n")

# 3. autocorrelation
# aim: no long term autocorr, 
draws_arr <- as.array(fit1)
bayesplot::mcmc_acf(draws_arr, pars = "alpha", lags = 30)


# 4. posterior‑predictive fit 
posterior <- rstan::extract(fit1)
y_rep <- posterior$y_rep
ppc_dens_overlay(y = stan_data_1$count, yrep = y_rep[1:200, ])


## marginal distribution check
ppc_dens_overlay(
  y    = stan_data_1$count,          # observed counts
  yrep = y_rep[sample(4000, 200), ] # 200 random draws; keep it legible
)

## For discrete data ➡ better: rootogram
ppc_rootogram(
  y    = stan_data_1$count,
  yrep = y_rep[sample(4000, 50), ]  # 50 draws is usually enough
)

ppc_stat(y = stan_data_1$count, yrep = y_rep, stat = "mean")



#7. credible interval
# too conservative?
y_pred_ci <- apply(y_rep, 2, quantile, probs = c(0.025, 0.975))
mean(stan_data_1$count >= y_pred_ci[1, ] & stan_data_1$count <= y_pred_ci[2, ])  # coverage rate



########  model1 forecasting check:
y_test = agg_df$Count[(N_train + 1):N_total]
y_test

fc_pois <- as_draws_matrix(rstan::extract(fit1, pars = "y_fore")$y_fore)

rmse <- function(mat, obs) {
  sqrt( mean( (colMeans(mat) - obs)^2 ) )
}
mae  <- function(mat, obs) {
  mean( abs(colMeans(mat) - obs) )
}

metrics <- tibble(
  model = c("Poisson‑SSM", "NegBin‑Region"),
  RMSE  = c(rmse(fc_pois, y_test), 0), #stub for model2's prediction
  MAE   = c(mae(fc_pois, y_test), 0)
)


crps_pois <- mean(crps_sample(y_test, t(fc_pois)))
crps_nb   <- 1
logs_pois <- mean(logs_sample(y_test, t(fc_pois)))
logs_nb   <- 1

metrics <- metrics |>
  mutate(
    CRPS = c(crps_pois, crps_nb),
    LogS = c(logs_pois, logs_nb)
  )

print(metrics)

pi95_pois <- mean(
  y_test >= apply(fc_pois, 2, quantile, 0.025) &
    y_test <= apply(fc_pois, 2, quantile, 0.975)
)


####### END OF model 1: Poisson model ######################
###########################################################################

###########################################################################
####### SECTION: model 2: per region, NegBinom model ######################

####### END OF model 2: per region, NegBinom model ######################
###########################################################################


###########################################################################
####### SECTION: model comparison ######################

####### END OF SECTION: model comparison ######################
###########################################################################


agg_df
