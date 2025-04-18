library(dplyr)
library(lubridate)
library(rstan)
library(tidyr)
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

stan_data <- list(
  N = nrow(agg_df),
  R = length(unique(df$region)),
  region = as.integer(factor(agg_df$region)),
  count = agg_df$Count,
  time = as.vector(agg_df$TimeIndex),
  depth = as.vector(agg_df$Depth),
  lat = as.vector(agg_df$Lat),
  lon = as.vector(agg_df$Lon),
  lag_mag = as.vector(agg_df$LagMag),
  nst = as.vector(agg_df$Nst),
  rms = as.vector(agg_df$RMS),
  clo = as.vector(agg_df$Clo)
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

fit <- stan(
  file = "RScripts/NegBin.stan",
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  init = "random"
)

posterior <- rstan::extract(fit)
y_rep <- posterior$y_rep

y_rep

hist(agg_df$Count, breaks = 30, col = rgb(0, 0, 1, 0.5),
     main = "Posterior Predictive Check", xlab = "Earthquake Counts")
lines(density(y_rep[1,]), col = "red")
