getwd()

library(ggplot2)
library(dplyr)
library(lubridate)

df <- read.csv("../cleaned_earthquake_data.csv")

# Combine Date and Time into POSIXct datetime
df$datetime <- as.POSIXct(paste(df$Date, df$Time), format="%Y/%m/%d %H:%M:%OS", tz="UTC")

# Drop rows with invalid datetime
df <- df[!is.na(df$datetime), ]

# Summary
summary(df)
colSums(is.na(df))

# latitude vs longitude: San Andreas Fault confirmed :)
ggplot(df, aes(x = Lon, y = Lat)) +
  geom_point(alpha = 0.3, color = "darkred") +
  labs(
    title = "Earthquake Epicenters in California",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

# Magnitude histogram
# this could give prior intuition?
ggplot(df, aes(x = Mag)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Earthquake Magnitude Distribution", x = "Magnitude", y = "Count")

# Depth histogram
# depth can influence magnitude: https://www.nature.com/articles/365045a0
ggplot(df, aes(x = Depth)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "black") +
  labs(title = "Earthquake Depth Distribution", x = "Depth (km)", y = "Count")


# Magnitudes type histogram:
# hierarchical model?
ggplot(df, aes(x = Magt)) +
  geom_bar(fill = "plum") +
  labs(title = "Magnitude Type Frequency", x = "Magnitude Type", y = "Count")

# Depth vs Magnitude for each magnitude type
ggplot(df, aes(x = Depth, y = Mag, color = Magt)) +
  geom_point(alpha = 0.5) +
  labs(
    title = "Depth vs Magnitude by Magnitude Type",
    x = "Depth (km)",
    y = "Magnitude",
    color = "Magnitude Type"
  ) +
  theme_minimal()

ggplot(df, aes(x = Depth, y = Mag)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  facet_wrap(~ Magt) +
  labs(
    title = "Depth vs Magnitude by Magnitude Type",
    x = "Depth (km)",
    y = "Magnitude"
  ) +
  theme_minimal()

# Nst vs RMS
# these two can indicate reliability of measurements, which can be used 
# as a variable to model the variance of magnitudes?
ggplot(df, aes(x = Nst, y = RMS)) +
  geom_point(alpha = 0.4, color = "purple") +
  coord_cartesian(ylim = c(0, 1)) +  # Zoom in on y-axis
  labs(
    title = "RMS Error vs Number of Seismic Stations (Zoomed)",
    x = "Number of Seismic Stations (Nst)",
    y = "RMS Error (zoomed to 0–1)"
  ) +
  theme_minimal()

# Earthquake counts per day, average magnitude per day, max magnitude per day
df_daily_summary <- df %>%
  mutate(date_only = as.Date(datetime)) %>%
  group_by(date_only) %>%
  summarise(
    daily_count = n(),
    avg_mag = mean(Mag, na.rm = TRUE),
    max_mag = max(Mag, na.rm = TRUE)
  )

p1 <- ggplot(df_daily_summary, aes(x = date_only, y = daily_count)) +
  geom_line(color = "steelblue") +
  labs(title = "Daily Earthquake Count", x = "Date", y = "Number of Earthquakes") +
  theme_minimal()

p2 <- ggplot(df_daily_summary, aes(x = date_only, y = avg_mag)) +
  geom_line(color = "darkgreen") +
  labs(title = "Daily Average Magnitude", x = "Date", y = "Average Magnitude") +
  theme_minimal()

p3 <- ggplot(df_daily_summary, aes(x = date_only, y = max_mag)) +
  geom_line(color = "firebrick") +
  labs(title = "Daily Maximum Magnitude", x = "Date", y = "Max Magnitude") +
  theme_minimal()

print(p1)
print(p2)
print(p3)

# per month...
df_monthly_summary <- df %>%
  mutate(month = floor_date(datetime, unit = "month")) %>%
  group_by(month) %>%
  summarise(
    monthly_count = n(),
    avg_mag = mean(Mag, na.rm = TRUE),
    max_mag = max(Mag, na.rm = TRUE)
  )

p1 <- ggplot(df_monthly_summary, aes(x = month, y = monthly_count)) +
  geom_line(color = "steelblue") +
  labs(title = "Monthly Earthquake Count", x = "Month", y = "Number of Earthquakes") +
  theme_minimal()


p2 <- ggplot(df_monthly_summary, aes(x = month, y = avg_mag)) +
  geom_line(color = "darkgreen") +
  labs(title = "Monthly Average Magnitude", x = "Month", y = "Average Magnitude") +
  theme_minimal()

p3 <- ggplot(df_monthly_summary, aes(x = month, y = max_mag)) +
  geom_line(color = "firebrick") +
  labs(title = "Monthly Maximum Magnitude", x = "Month", y = "Max Magnitude") +
  theme_minimal()

print(p1)
print(p2)
print(p3)

# per year...
df_yearly_summary <- df %>%
  mutate(year = year(datetime)) %>%
  group_by(year) %>%
  summarise(
    yearly_count = n(),
    avg_mag = mean(Mag, na.rm = TRUE),
    max_mag = max(Mag, na.rm = TRUE)
  )

p1 <- ggplot(df_yearly_summary, aes(x = year, y = yearly_count)) +
  geom_line(color = "steelblue") +
  labs(title = "Yearly Earthquake Count", x = "Year", y = "Number of Earthquakes") +
  theme_minimal()

p2 <- ggplot(df_yearly_summary, aes(x = year, y = avg_mag)) +
  geom_line(color = "darkgreen") +
  labs(title = "Yearly Average Magnitude", x = "Year", y = "Average Magnitude") +
  theme_minimal()

p3 <- ggplot(df_yearly_summary, aes(x = year, y = max_mag)) +
  geom_line(color = "firebrick") +
  labs(title = "Yearly Maximum Magnitude", x = "Year", y = "Max Magnitude") +
  theme_minimal()

print(p1)
print(p2)
print(p3)

# we should do analysis on month or year basis





