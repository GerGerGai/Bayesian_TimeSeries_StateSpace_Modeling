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

stan_data <- list(
  N = nrow(agg_df),
  count = agg_df$Count,
  time = as.vector(agg_df$TimeIndex),
  depth = as.vector(agg_df$Depth),
  lag_mag = as.vector(agg_df$LagMag)
)


init_fn <- function() {
  R <- stan_data$R
  
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

fit <- stan(
  file = "RScripts/Poisson_simple.stan",
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 42,
  control = list(adapt_delta = 0.8, max_treedepth = 10),
  init = "random"
)

# check mixing/convergence:
print(fit)
summary(fit)$summary

#1. check fitness of the model: the result seems good, though not perfect
posterior <- rstan::extract(fit)
y_rep <- posterior$y_rep
ppc_dens_overlay(y = stan_data$count, yrep = y_rep[1:200, ])


#marginal distribution check
ppc_dens_overlay(
  y    = stan_data$count,          # observed counts
  yrep = y_rep[sample(4000, 200), ] # 200 random draws; keep it legible
)

## For discrete data ➡ better: rootogram
ppc_rootogram(
  y    = stan_data$count,
  yrep = y_rep[sample(4000, 50), ]  # 50 draws is usually enough
)



#ci
y_pred_ci <- apply(y_rep, 2, quantile, probs = c(0.05, 0.95))
mean(stan_data$count >= y_pred_ci[1, ] & stan_data$count <= y_pred_ci[2, ])  # coverage rate


