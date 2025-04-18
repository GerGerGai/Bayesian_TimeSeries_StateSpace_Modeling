library(dplyr)
library(lubridate)
library(rstan)
library(tidyr)

## Data Preparation
df <- read.csv("cleaned_earthquake_data.csv")
df$datetime <- as.POSIXct(paste(df$Date, df$Time), format="%Y/%m/%d %H:%M:%OS", tz="UTC") # combine dates and times

df <- df %>%
  filter(!is.na(datetime) & Mag >= 4) %>%
  mutate(
    month = floor_date(datetime, "month"), # convert to month level
    region = paste0("Lat", floor(Lat), "_Lon", floor(Lon)) # define grid regions of 1°x1°
  )

agg_df <- df %>%
  group_by(month, region) %>%
  summarise(
    Count = n(),
    Depth = mean(Depth, na.rm = TRUE),
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE),
    Mag = mean(Mag, na.rm = TRUE),
    Nst = mean(Nst, na.rm = TRUE),
    RMS = mean(RMS, na.rm = TRUE),
    Clo = mean(Clo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month, region) %>%
  group_by(region) %>%
  mutate(
    LagMag = lag(Mag), # previous month's avg magnitude
    TimeIndex = dense_rank(month) # numeric time feature for Stan
  ) %>%
  ungroup() %>%
  drop_na()


stan_data <- list(
  N = nrow(agg_df),
  R = length(unique(agg_df$region)),
  region = as.integer(factor(agg_df$region)), # convert regions into region ID
  y = agg_df$Count,
  depth = agg_df$Depth,
  lat = agg_df$Lat,
  lon = agg_df$Lon,
  time = agg_df$TimeIndex,
  lag_mag = agg_df$LagMag,
  nst = agg_df$Nst,
  rms = agg_df$RMS,
  clo = agg_df$Clo
)

init_fn <- function() {
  list(
    beta = rep(0, 6),         
    gamma = rep(0, 3),       
    sigma_u = 1,
    u_raw = rnorm(stan_data$R)
  )
}

fit <- stan(
  file = "RScripts/NegBin.stan",
  data = stan_data,
  iter = 2000, 
  chains = 4, 
  seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  init = init_fn
)

## Posterior Predictive Checks
posterior <- rstan::extract(fit)
y_rep <- posterior$y_rep

# Histogram overlay
hist(agg_df$Count, breaks=30, col=rgb(0,0,1,0.5), main="Posterior Predictive Check", xlab="Earthquake Counts")
lines(density(y_rep[1,]), col="red")


