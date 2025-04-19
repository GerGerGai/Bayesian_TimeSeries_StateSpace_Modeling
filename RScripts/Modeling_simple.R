library(dplyr)
library(lubridate)
library(rstan)
library(tidyr)
library(bayesplot)
#install.packages("V8")


## Data Preparation
df <- read.csv("cleaned_earthquake_data.csv")
df$datetime <- as.POSIXct(paste(df$Date, df$Time), format="%Y/%m/%d %H:%M:%OS", tz="UTC") # combine dates and times


df <- df %>%
  filter(!is.na(datetime) & Mag >= 3.5) 

df <- df %>%
  mutate(LagMag = lag(Mag)) %>%
  ungroup()



# Time index as continuous covariate (in months since first quake)
origin_time <- min(df$datetime)
df$TimeIndex <- as.numeric(difftime(df$datetime, origin_time, units = "days")) / 30




agg_df <- df %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarise(
    Count = n(),
    Depth = mean(Depth, na.rm = TRUE),
    LagMag = mean(LagMag, na.rm = TRUE),
    TimeIndex = mean(TimeIndex, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  drop_na()

agg_df <- agg_df %>%
  mutate(across(c(Depth, LagMag, TimeIndex), scale)) %>%
  mutate(across(everything(), ~ifelse(is.finite(.), ., 1e-3)))


###########################################################################
####### SECTION: model 1: Poisson model ######################

stan_data_1 <- list(
  N = nrow(agg_df),
  count = agg_df$Count,
  time = as.vector(agg_df$TimeIndex),
  depth = as.vector(agg_df$Depth),
  lag_mag = as.vector(agg_df$LagMag)
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
y_pred_ci <- apply(y_rep, 2, quantile, probs = c(0.05, 0.95))
mean(stan_data_1$count >= y_pred_ci[1, ] & stan_data_1$count <= y_pred_ci[2, ])  # coverage rate



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
