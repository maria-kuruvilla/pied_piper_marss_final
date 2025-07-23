# Goal - to create a figure of raw data for the manuscript with last 5 years of 
# counts from all three rivers

# Load libraries

library(ggplot2)
library(here)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(zoo)

# Load data


data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"

data_day_night <- rbind(data_day,data_night)

data_day_night_subset_long <- data_day_night %>% 
  filter(daytime_category == "Night") %>% 
  dplyr::select(year, doy, chinook0_hatchery_perhour, chinook0_wild_perhour, 
         chinook0_hatchery_perhour_interpolate) %>% 
  pivot_longer(cols = c(chinook0_wild_perhour, chinook0_hatchery_perhour, 
                        chinook0_hatchery_perhour_interpolate), 
               names_to = c("origin" , "interpolate"),
               names_pattern = "chinook0_(.*)_perhour(.*)",
               values_to = "fish_perhour") %>% 
  mutate(origin_interpolate = ifelse(interpolate=="", origin,"hatchery\ninterpolated"))


w1_new <- ggplot(data_day_night_subset_long %>% 
                   filter(year >=2016) %>% 
                   filter(doy > 130, doy < 200)) +
  geom_line(aes(doy, y = fish_perhour, color = origin_interpolate, 
                linetype = origin_interpolate),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("wild" = "salmon","hatchery" = "cadetblue",
                                "hatchery\ninterpolated" = "cadetblue"))+
  scale_linetype_manual(values = c("wild" = "solid","hatchery" = "solid",
                                   "hatchery\ninterpolated" = "twodash"))+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 140,180))+
  theme_classic()+
  labs(x = "Day of Year", y = "Number of Chinook per hour")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(l = 10)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))

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
  scale_y_continuous(n.breaks = 4) +
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

data_subset_long <- data %>% 
  # filter(daytime_category == "Night") %>% 
  dplyr::select(year, doy, chinook0_hatchery_perhour_night, chinook0_wild_perhour_night) %>%
  mutate(chinook0_hatchery_perhour_night_interpolate = chinook0_hatchery_perhour_night) %>% 
  mutate(chinook0_hatchery_perhour_night = ifelse(is.na(chinook0_wild_perhour_night), NA, chinook0_hatchery_perhour_night_interpolate)) %>% 
  group_by(year,doy) %>% 
  # View() %>% 
  pivot_longer(cols = c(chinook0_wild_perhour_night, chinook0_hatchery_perhour_night,
                        chinook0_hatchery_perhour_night_interpolate), 
               names_to = c("origin" , "interpolate"),
               names_pattern = "chinook0_(.*)_perhour_night(.*)",
               values_to = "fish_perhour") %>% 
  mutate(origin_interpolate = ifelse(interpolate=="", origin,"hatchery\ninterpolated"))

w2_new <- ggplot(data_subset_long %>% 
                   filter(year >=2016, year <=2020) %>% 
                   filter(doy > 130, doy < 200)) +
  geom_line(aes(doy, y = fish_perhour, color = origin_interpolate, 
                linetype = origin_interpolate),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("wild" = "salmon","hatchery" = "cadetblue",
                                "hatchery\ninterpolated" = "cadetblue"))+
  scale_linetype_manual(values = c("wild" = "solid","hatchery" = "solid",
                                   "hatchery\ninterpolated" = "twodash"))+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 140, 180, 220))+
  scale_y_continuous(n.breaks = 4) +
  theme_classic()+
  labs(x = "Day of Year", y = "Number of Chinook per hour")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(l = 10)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))

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

data_day_night_skagit_subset_long <- data_day_night %>% 
  filter(daytime_category == "night") %>% 
  dplyr::select(year, doy, chinook0_hatchery_perhour, chinook0_wild_perhour, 
         chinook0_hatchery_perhour_inp) %>% 
  pivot_longer(cols = c(chinook0_wild_perhour, chinook0_hatchery_perhour, 
                        chinook0_hatchery_perhour_inp), 
               names_to = c("origin" , "interpolate"),
               names_pattern = "chinook0_(.*)_perhour(.*)",
               values_to = "fish_perhour") %>% 
  mutate(origin_interpolate = ifelse(interpolate=="", origin,"hatchery\ninterpolated"))

w3_new <- ggplot(data_day_night_skagit_subset_long %>% 
                   filter(year >=2016, year<=2020) %>% 
                   filter(doy > 130, doy < 200)) +
  geom_line(aes(doy, y = fish_perhour, color = origin_interpolate, 
                linetype = origin_interpolate),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("wild" = "salmon","hatchery" = "cadetblue",
                                "hatchery\ninterpolated" = "cadetblue"))+
  scale_linetype_manual(values = c("wild" = "solid","hatchery" = "solid",
                                   "hatchery\ninterpolated" = "twodash"))+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 140,180))+
  theme_classic()+
  labs(x = "Day of Year", y = "Number of Chinook per hour")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(l = 10)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))


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

(w1_new/w2_new/w3_new 
  + plot_layout(guides = "collect", axis_titles = "collect")+ plot_annotation(tag_levels = 'a'))

ggsave(here("visualizations","output","chinook_data_raw_five_years_new.png"), width = 10, height = 8, 
       units = "in", dpi = 300)


# coho

#interpolate the hatchery covariate values
data_day$coho1_hatchery_perhour_interpolate <- na.approx(data_day$coho1_hatchery_perhour, na.rm = FALSE)
data_night$coho1_hatchery_perhour_interpolate <- na.approx(data_night$coho1_hatchery_perhour, na.rm = FALSE)


data_day_night <- rbind(data_day,data_night)

data_day_night_coho_subset_long <- data_day_night %>% 
  filter(daytime_category == "Night") %>% 
  select(year, doy, coho1_hatchery_perhour, coho1_wild_perhour, 
         coho1_hatchery_perhour_interpolate) %>% 
  pivot_longer(cols = c(coho1_wild_perhour, coho1_hatchery_perhour, 
                        coho1_hatchery_perhour_interpolate), 
               names_to = c("origin" , "interpolate"),
               names_pattern = "coho1_(.*)_perhour(.*)",
               values_to = "fish_perhour") %>% 
  mutate(origin_interpolate = ifelse(interpolate=="", origin,"hatchery\ninterpolated"))


w1_coho_new <- ggplot(data_day_night_coho_subset_long %>% 
                   filter(year >=2016) %>% 
                   filter(doy > 100, doy < 160)) +
  geom_line(aes(doy, y = fish_perhour, color = origin_interpolate, 
                linetype = origin_interpolate),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("wild" = "salmon","hatchery" = "cadetblue",
                                "hatchery\ninterpolated" = "cadetblue"))+
  scale_linetype_manual(values = c("wild" = "solid","hatchery" = "solid",
                                   "hatchery\ninterpolated" = "twodash"))+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100,150))+
  theme_classic()+
  labs(x = "Day of Year", y = "Number of coho per hour")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(l = 10)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))

# puyallup

data_coho_subset_long <- data %>% 
  # filter(daytime_category == "Night") %>% 
  select(year, doy, coho1_hatchery_perhour_night, coho1_wild_perhour_night) %>%
  mutate(coho1_hatchery_perhour_night_interpolate = coho1_hatchery_perhour_night) %>% 
  mutate(coho1_hatchery_perhour_night = ifelse(is.na(coho1_wild_perhour_night), NA, coho1_hatchery_perhour_night_interpolate)) %>% 
  group_by(year,doy) %>% 
  # View() %>% 
  pivot_longer(cols = c(coho1_wild_perhour_night, coho1_hatchery_perhour_night,
                        coho1_hatchery_perhour_night_interpolate), 
               names_to = c("origin" , "interpolate"),
               names_pattern = "coho1_(.*)_perhour_night(.*)",
               values_to = "fish_perhour") %>% 
  mutate(origin_interpolate = ifelse(interpolate=="", origin,"hatchery\ninterpolated"))

w2_coho_new <- ggplot(data_coho_subset_long %>% 
                   filter(year >=2016, year <=2020) %>% 
                   filter(doy > 100, doy < 160)) +
  geom_line(aes(doy, y = fish_perhour, color = origin_interpolate, 
                linetype = origin_interpolate),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("wild" = "salmon","hatchery" = "cadetblue",
                                "hatchery\ninterpolated" = "cadetblue"))+
  scale_linetype_manual(values = c("wild" = "solid","hatchery" = "solid",
                                   "hatchery\ninterpolated" = "twodash"))+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 150))+
  scale_y_continuous(n.breaks = 4) +
  theme_classic()+
  labs(x = "Day of Year", y = "Number of coho per hour")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(l = 10)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))



skagit_data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))


skagit_data_day_night$coho1_wild_perhour <- skagit_data_day_night$coho1_wild_num/skagit_data_day_night$In
skagit_data_day_night$coho1_hatchery_perhour <- skagit_data_day_night$coho1_hatchery_num/skagit_data_day_night$In


skagit_data_day_night$coho1_hatchery_perhour_inp <- na.approx(skagit_data_day_night$coho1_hatchery_perhour, 
                                                          na.rm = FALSE)

data_day_night_skagit_subset_long <- skagit_data_day_night %>% 
  filter(daytime_category == "night", trap == "scoop") %>% 
  dplyr::select(year, doy, coho1_hatchery_perhour, coho1_wild_perhour, 
         coho1_hatchery_perhour_inp) %>% 
  pivot_longer(cols = c(coho1_wild_perhour, coho1_hatchery_perhour, 
                        coho1_hatchery_perhour_inp), 
               names_to = c("origin" , "interpolate"),
               names_pattern = "coho1_(.*)_perhour(.*)",
               values_to = "fish_perhour") %>% 
  mutate(origin_interpolate = ifelse(interpolate=="", origin,"hatchery\ninterpolated"))

w3_coho_new <- ggplot(data_day_night_skagit_subset_long %>% 
                   filter(year >=2016, year<=2020) %>% 
                   filter(doy > 100, doy < 160)) +
  geom_line(aes(doy, y = fish_perhour, color = origin_interpolate, 
                linetype = origin_interpolate),
            linewidth = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("wild" = "salmon","hatchery" = "cadetblue",
                                "hatchery\ninterpolated" = "cadetblue"))+
  scale_linetype_manual(values = c("wild" = "solid","hatchery" = "solid",
                                   "hatchery\ninterpolated" = "twodash"))+
  facet_wrap(~year, ncol = 5, scales = "free")+
  scale_x_continuous(breaks = c(100, 150))+
  theme_classic()+
  labs(x = "Day of Year", y = "Number of coho per hour")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(l = 10)),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 14))

w3_coho_new

(w1_coho_new/w2_coho_new/w3_coho_new + plot_layout(guides = "collect", axis_titles = "collect")+ plot_annotation(tag_levels = 'a'))

ggsave(here("visualizations","output","coho_data_raw_five_years_new.png"), width = 10, height = 8, 
       units = "in", dpi = 300)

(w1_coho_new/w2_coho_new/w3_coho_new + plot_annotation(tag_levels = 'a')+
    plot_layout(guides = "collect", axis_titles = "collect"))

ggsave(here("visualizations","output","coho_data_raw_five_years_new.png"), width = 10, height = 8, 
       units = "in", dpi = 300)








