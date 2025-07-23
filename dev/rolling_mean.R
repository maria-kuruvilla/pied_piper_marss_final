# calculate rolling mean of flow, temp, etc and then subtract it from daily data to get anomaly

# Load required libraries

library(tidyverse)
library(ggplot2)
library(here)
library(zoo)

# Load data

data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

#interpolate the hatchery covariate values
data_day$coho1_hatchery_perhour_interpolate <- na.approx(data_day$coho1_hatchery_perhour, na.rm = FALSE)
data_night$coho1_hatchery_perhour_interpolate <- na.approx(data_night$coho1_hatchery_perhour, na.rm = FALSE)


data_day_night <- rbind(data_day,data_night)

#calculate the cumulative sum of chinook0_wild_perhour

data_day_night <- data_day_night %>%
  filter(doy > 120, doy < 250) %>%
  group_by(year, daytime_category) %>% 
  mutate(chinook0_wild_perhour_cumsum = cumsum(ifelse(is.na(chinook0_wild_perhour),
                                                      0,
                                                      chinook0_wild_perhour))) %>% 
  mutate(chinook0_proportion = chinook0_wild_perhour_cumsum/
           sum(chinook0_wild_perhour,na.rm = TRUE)) %>% 
  mutate(temp_rolling_mean = rollmean(temp, k = 32, fill = NA, align = "right"),
         flow_rolling_mean = rollmean(flow, k = 32, fill = NA, align = "right"))


ggplot(data_day_night %>% filter(daytime_category == "Day") %>% 
         mutate(temp_rolling_mean = rollmean(temp, k = 32, fill = NA, align = "right"),
                flow_rolling_mean = rollmean(flow, k = 32, fill = NA, align = "right"),
                temp_anomaly = temp - temp_rolling_mean,
                flow_anomaly = flow - flow_rolling_mean)) +
  geom_line(aes(x = doy, y = temp, color = daytime_category)) +
  geom_line(aes(x = doy, y = temp_rolling_mean), color = "black") +
  geom_line(aes(x = doy, y = temp_anomaly), color = "salmon")+
  facet_wrap(~year) 

# plot the difference in chinook0_hatchery_perhour for day and night

data_day_night %>% 
  group_by(year,daytime_category) %>% 
  mutate(chinook0_hatchery_perhour_diff = c(NA,diff(chinook0_hatchery_perhour_interpolate,1))) %>% 
  mutate(coho1_hatchery_perhour_diff = c(NA,diff(coho1_hatchery_perhour_interpolate,1))) %>%
  # mutate(chinook0_hatchery_perhour_diff_interpolate = na.approx(chinook0_hatchery_perhour_diff, na.rm = FALSE)) %>% 
  filter(doy > 100, doy < 200, daytime_category == "Night", year== 2016) %>% 
  ggplot()+
  geom_line(aes(x = doy, y = coho1_hatchery_perhour_diff), color = "cyan", alpha = 0.5, size = 1.2)+
  # geom_line(aes(x = doy, y = chinook0_hatchery_perhour_diff_interpolate), color = "cyan", alpha = 0.5, size = 1.2, linetype = "dotted")+
  geom_line(aes(x = doy , y = coho1_hatchery_perhour), color = "salmon", alpha = 0.5, size = 1.2) +
  geom_line(aes(x = doy , y = coho1_hatchery_perhour_interpolate), color = "salmon", alpha = 0.5, size = 1.2, linetype = "dotted") +
  facet_wrap(~year, scales = "free")+
  theme_classic()

#calculate the rolling mean and the hatchery_difference variables and save the csv

data_day_night_new <- data_day_night %>% 
  group_by(year,daytime_category) %>% 
  mutate(chinook0_hatchery_perhour_diff = c(NA,diff(chinook0_hatchery_perhour_interpolate,1))) %>% 
  mutate(coho1_hatchery_perhour_diff = c(NA,diff(coho1_hatchery_perhour_interpolate,1))) %>%
  mutate(temp_rolling_mean = zoo::rollapply(temp, width = 32, FUN = mean, na.rm = TRUE, fill = NA, align = "right"),
         flow_rolling_mean = zoo::rollapply(flow, width = 32,  FUN = mean, na.rm = TRUE, fill = NA, align = "right"),
         temp_anomaly = temp - temp_rolling_mean,
         flow_anomaly = flow - flow_rolling_mean)


# plot temp_anomaly for 2016 from data_day_night_new

data_day_night_new %>% 
  filter(year == 2016, doy > 100, doy < 200) %>% 
  ggplot()+
  geom_line(aes(x = doy, y = temp), color = "black")+
  geom_line(aes(x = doy, y = temp_anomaly, color = daytime_category))+
  facet_wrap(~year)

#save csv

write.csv(data_day_night_new,here("data","dungeness_anomaly.csv"))
