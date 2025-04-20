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
  filter(!is.na(datetime) & Mag >= 4) 

region_size <- 2
df$region_lat <- floor(df$Lat / region_size)
df$region_lon <- floor(df$Lon / region_size)
df$region <- paste0("Lat", df$region_lat, "_Lon", df$region_lon)
df$region_idx <- as.numeric(as.factor(df$region))

# Time index as continuous covariate (in months since first quake)
origin_time <- min(df$datetime)
df$TimeIndex <- as.numeric(difftime(df$datetime, origin_time, units = "days")) / 30

# Add sine/cosine seasonal terms for monthly seasonality
df$SeasonSin <- sin(2 * pi * df$TimeIndex / 12)
df$SeasonCos <- cos(2 * pi * df$TimeIndex / 12)

df <- df %>%
  arrange(region, TimeIndex) %>%
  group_by(region) %>%
  mutate(LagMag = lag(Mag)) %>%
  ungroup()


agg_df <- df %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(region, month) %>%
  summarise(
    Count = n(),
    Depth = mean(Depth, na.rm = TRUE),
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE),
    LagMag = mean(LagMag, na.rm = TRUE),
    TimeIndex = mean(TimeIndex, na.rm = TRUE),
    SeasonSin = mean(SeasonSin, na.rm = TRUE),
    SeasonCos = mean(SeasonCos, na.rm = TRUE),
    Nst = mean(Nst, na.rm = TRUE),
    RMS = mean(RMS, na.rm = TRUE),
    Clo = mean(Clo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  drop_na()

agg_df <- agg_df %>%
  mutate(across(c(Depth, Lat, Lon, LagMag, TimeIndex, SeasonSin, SeasonCos, Nst, RMS, Clo), scale)) %>%
  mutate(across(everything(), ~ifelse(is.finite(.), ., 1e-3)))

set.seed(42)
n <- nrow(agg_df)
train_idx <- sample(1:n, size = floor(0.8 * n), replace = FALSE)
agg_train <- agg_df[train_idx, ]
agg_test <- agg_df[-train_idx, ]

stan_data <- list(
  N = nrow(agg_train),
  R = length(unique(df$region)),
  region = as.integer(factor(agg_train$region)),
  count = agg_train$Count,
  time = as.vector(agg_train$TimeIndex),
  depth = as.vector(agg_train$Depth),
  lat = as.vector(agg_train$Lat),
  lon = as.vector(agg_train$Lon),
  lag_mag = as.vector(agg_train$LagMag),
  nst = as.vector(agg_train$Nst),
  rms = as.vector(agg_train$RMS),
  clo = as.vector(agg_train$Clo)
)

init_fn <- function() {
  R <- stan_data$R
  
  list(
    alpha = rnorm(R, 0, 0.1),
    beta_time = rnorm(R, 0, 0.1),
    beta_sin = rnorm(R, 0, 0.1),
    beta_cos = rnorm(R, 0, 0.1),
    
    theta_depth = 0,
    theta_lat = 0,
    theta_lon = 0,
    theta_lagmag = 0,
    
    gamma = rep(0, 3),
    
    mu_alpha = 0,
    sigma_alpha = 1,
    
    mu_beta_time = 0,
    sigma_beta_time = 1,
    
    mu_beta_sin = 0,
    sigma_beta_sin = 1,
    
    mu_beta_cos = 0,
    sigma_beta_cos = 1
  )
}

options(mc.cores = 4)

fit <- stan(
  file = "RScripts/NegBin.stan",
  data = stan_data,
  iter = 1000,
  chains = 4,
  seed = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  init = "random"
)

# save the model to disk first:
saveRDS(fit, file = "trained_model_tail.rds")
# to read: fit <- readRDS("trained_model.rds")

# check mixing/convergence:
print(fit)
summary(fit)$summary

# check fitness of the model: the result seems good, though not perfect
posterior <- rstan::extract(fit)
y_rep <- posterior$y_rep
ppc_dens_overlay(y = stan_data$count, yrep = y_rep[1:100, ])

y_pred_ci <- apply(y_rep, 2, quantile, probs = c(0.05, 0.95))
mean(stan_data$count >= y_pred_ci[1, ] & stan_data$count <= y_pred_ci[2, ])  

y_pred_mean <- colMeans(y_rep)
plot(y_pred_mean, stan_data$count,
     xlab = "Predicted mean", ylab = "Observed count", main = "Fit scatterplot")
abline(0, 1, col = "red")

resid <- stan_data$count - y_pred_mean
plot(resid, main = "Residuals", ylab = "Observed - Predicted")
abline(h = 0, col = "red", lty = 2)

region_test <- as.integer(factor(agg_test$region))  
time_test   <- agg_test$TimeIndex
depth_test  <- agg_test$Depth
lat_test    <- agg_test$Lat
lon_test    <- agg_test$Lon
lagmag_test <- agg_test$LagMag
nst_test    <- agg_test$Nst
rms_test    <- agg_test$RMS
clo_test    <- agg_test$Clo

n_test <- nrow(agg_test)
n_draws <- length(posterior$theta_depth)
y_rep_test <- matrix(NA, nrow = n_draws, ncol = n_test)

for (d in 1:n_draws) {
  eta <- posterior$alpha[d, region_test] +
    posterior$beta_time[d, region_test] * time_test +
    posterior$beta_sin[d, region_test] * sin(2 * pi * time_test / 12) +
    posterior$beta_cos[d, region_test] * cos(2 * pi * time_test / 12) +
    posterior$theta_depth[d] * depth_test +
    posterior$theta_lat[d]   * lat_test +
    posterior$theta_lon[d]   * lon_test +
    posterior$theta_lagmag[d] * lagmag_test
  
  mu <- exp(eta)
  
  phi <- exp(
    posterior$gamma[d, 1] * nst_test +
      posterior$gamma[d, 2] * rms_test +
      posterior$gamma[d, 3] * clo_test
  )
  
  y_rep_test[d, ] <- rnbinom(n_test, size = phi, mu = mu)
}

y_pred_mean <- colMeans(y_rep_test)

plot(y_pred_mean, agg_test$Count,
     xlab = "Predicted Mean", ylab = "Observed Count",
     main = "Test Set Prediction")
abline(0, 1, col = "red")

rmse <- sqrt(mean((y_pred_mean - agg_test$Count)^2))

y_lower <- apply(y_rep_test, 2, quantile, probs = 0.025)
y_upper <- apply(y_rep_test, 2, quantile, probs = 0.975)
coverage <- mean(agg_test$Count >= y_lower & agg_test$Count <= y_upper)

cat("Test RMSE:", round(rmse, 2), "\n")
cat("95% Coverage:", round(coverage * 100, 2), "%\n")



