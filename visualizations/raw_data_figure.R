# Goal - to create a figure of raw data for the manuscript with last 5 years of 
# counts from all three rivers

# Load libraries

library(ggplot2)
library(here)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Load data


data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)

w1 <- ggplot() +
  geom_line(data = data_day_night %>% 
              arrange(doy ,year) %>%
              filter(daytime_category == 'Night',year>= 2016,
                     doy >= 100 & doy <= 200),
            aes(x = doy, y = scale(chinook0_wild_perhour), color = "wild"), 
            linewidth = 1.2, alpha = 0.5)+
  geom_line(data = data_day_night %>% 
              arrange(doy ,year) %>%
              filter(daytime_category == 'Night',year >= 2016 ,
                     doy >= 100 & doy <= 200),
            aes(x = doy, y = scale(chinook0_hatchery_perhour, center = FALSE), color = "hatchery"), 
            linewidth = 1.2, alpha = 0.5)+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 140,180))+
  scale_color_manual(name = "",labels = c("wild","hatchery"),values = c("wild" = "salmon","hatchery" = "cadetblue"))+
  # scale_alpha_manual(name = "", labels = c("hatchery", "wild"),values = c("wild" = 0.5,"hatchery" = 0.8))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled number of\n Chinook per hour", color = "Origin")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))



# puyallup

data <- read.csv(here("puyallup", "data","puyallup_final.csv"),
                 na.strings = c("NA",""))

w2 <- ggplot() +
  geom_line(data = data %>% 
              arrange(doy ,year) %>%
              filter(year>2016,
                     doy >= 100 & doy <= 200),
            aes(x = doy, y = scale(chinook0_wild_perhour_night), color = "wild"), 
            linewidth = 1.2, alpha = 0.5)+
  geom_line(data = data %>% 
              arrange(doy ,year) %>%
              filter(year>2016,
                     doy >= 100 & doy <= 200),
            aes(x = doy, y = scale(chinook0_hatchery_perhour_night, center = FALSE), color = "hatchery"), 
            linewidth = 1.2, alpha = 0.5)+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 140,180))+
  scale_color_manual(name = "",labels = c("wild","hatchery"),values = c("wild" = "salmon","hatchery" = "cadetblue"))+
  # scale_alpha_manual(name = "", labels = c("hatchery", "wild"),values = c("wild" = 0.5,"hatchery" = 0.8))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled number of\n Chinook per hour", color = "Origin")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))


#skagit


data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))


data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In


data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                          na.rm = FALSE)

w3 <- ggplot()+
  geom_line(data = data_day_night %>% 
              arrange(doy ,year) %>%
              filter(year>2017, daytime_category == 'night', trap == 'screw',
                     doy >= 100 & doy <= 200),
            aes(x = doy, y = scale(chinook0_wild_perhour), color = "wild"), 
            linewidth = 1.2, alpha = 0.5)+
  geom_line(data = data_day_night %>% 
              arrange(doy ,year) %>%
              filter(year>2017,daytime_category == 'night', trap == 'screw',
                     doy >= 100 & doy <= 200),
            aes(x = doy, y = scale(chinook0_hatchery_perhour, center = FALSE), color = "hatchery"), 
            linewidth = 1.2, alpha = 0.5)+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 140,180))+
  scale_color_manual(name = "",labels = c("wild","hatchery"),values = c("wild" = "salmon","hatchery" = "cadetblue"))+
  # scale_alpha_manual(name = "", labels = c("hatchery", "wild"),values = c("wild" = 0.5,"hatchery" = 0.8))+
  theme_classic()+
  labs(x = "Day of Year", y = "Scaled number of\n Chinook per hour", color = "Origin")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))

w1/w2/w3

ggsave(here("visualizations","output","chinook_data_raw_five_years.png"), width = 10, height = 10, 
       units = "in", dpi = 300)

