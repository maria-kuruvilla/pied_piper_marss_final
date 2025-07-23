


library(MARSS)
library(ggplot2)
library(here)
library(tidyverse)
library(ggpubr)
library(patchwork)



nyears = num_years*2
c = NULL
name_individual = NULL
for(kk in c(1,3,6)){
  c = rbind(c,covariates_coho1_skagit_night[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_coho1_skagit_night)[1+(kk-1)*nyears]
  name_individual = paste(name_individual, substr(name_long,1,nchar(name_long)-15))
  
}

print(name_individual)
fit.model = c(list(c= c), mod_list(nyears,3,1))
fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))


predict_x_coho1_skagit_w_hatchery<- predict(fit, type = "ytT", 
                                            interval = "confidence")


#predict ytT with new data, t will be 10, c will be rep(0,10), rep(0,10), rep(0,10), rep(0,10), rep(0,10)

c_new_data = matrix(0, nrow = 24*3, ncol = 10)
y_newdata = subset_coho_summer_perhour_night[,20:29]
timesteps <- dim(y_newdata)[2]
predict_no_hatchery <- predict(skagit_coho_best_model, newdata = list(y = y_newdata, t = 1:timesteps, c = c_new_data), type = "ytT", 
                               interval = "confidence")

predict_w_data <- predict(fit, newdata = list(y = y_newdata, 
                                              t = 1:timesteps, 
                                              c = c[,20:29]), 
                          type = "ytT", 
                          interval = "confidence")

plot(predict_no_hatchery)
plot(predict_w_data)
plot(predict(fit))



#predict x

predict_x_coho1_skagit_w_hatchery <- predict(fit, type = "xtT", 
                                             interval = "confidence")

predict_x_coho1_skagit_w_hatchery$pred$trap <- substring(unique(predict_x_coho1_skagit_w_hatchery$pred$.rownames), 14, 
                                                         length(unique(predict_x_coho1_skagit_w_hatchery$pred$.rownames)))

predict_x_coho1_skagit_w_hatchery$pred$year <- substring(unique(predict_x_coho1_skagit_w_hatchery$pred$.rownames), 3,6)

predict_x_coho1_skagit_w_hatchery$pred$daynight_category <- substring(unique(predict_x_coho1_skagit_w_hatchery$pred$.rownames), 8,12)

skagit_covariates_coho1 <- as.data.frame(t(covariates_coho1_skagit_night))

glimpse(skagit_covariates_coho1)
skagit_covariates_coho1$doy <- as.numeric(rownames(skagit_covariates_coho1))
skagit_covariates_coho1_long <- skagit_covariates_coho1 %>% 
  pivot_longer(cols = -c("doy"), names_to = c(".value","year", "daynight_category","trap"),
               names_pattern = "(.*)_(.{4})_(.{5})_(.{5})")
head(skagit_covariates_coho1_long)
#make column for doy


#merge skagit_covariates_coho1_long$coho1_hatchery_perhour_inp with predict_coho1_skagit$pred

predict_x_coho1_skagit_w_hatchery$pred2 <- predict_x_coho1_skagit_w_hatchery$pred %>% 
  mutate(doy = t+100) %>% 
  left_join(skagit_covariates_coho1_long, by = c("year", "daynight_category", "trap","doy"))




ggplot(data = predict_x_coho1_skagit_w_hatchery$pred %>% 
         filter(year == 2016 & daynight_category == "night" & trap == 'scoop'))+
  
  geom_line(aes(x = t+100, y = estimate, color = "predicted"))+
  geom_point(aes(x = t+100, y = .x, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = skagit_covariates_coho1, aes(x = doy, y = percent_hatchery_diff_2016_night_scoop), 
            color = "cadetblue", linewidth = 1.2, alpha = 0.7) +
  scale_color_manual(name = "Legend", values = c("predicted" = "salmon", "observed" = "salmon"))+
  labs(x = "Time", y = "log(scaled(coho salmon per hour) + 1)", title = "")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.2, linetype="solid"
    ),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )

scaled_data * attr(scaled_data, 'scaled:scale')

#predict without hatchery, use only first 2 covariates (flow anomaly and season)
c_new_data = rbind(c[1:48,],matrix(0, nrow = 24, ncol = dim(c)[2]))
y_newdata = subset_coho_summer_perhour_night

predict_x_coho1_skagit_wo_hatchery <- predict(fit, newdata = list(y = y_newdata, c = c_new_data), type = "xtT", 
                                              interval = "confidence")

predict_x_coho1_skagit_wo_hatchery$pred$trap <- substring(unique(predict_x_coho1_skagit_wo_hatchery$pred$.rownames), 14, 
                                                          length(unique(predict_x_coho1_skagit_wo_hatchery$pred$.rownames)))

predict_x_coho1_skagit_wo_hatchery$pred$year <- substring(unique(predict_x_coho1_skagit_wo_hatchery$pred$.rownames), 3,6)

predict_x_coho1_skagit_wo_hatchery$pred$daynight_category <- substring(unique(predict_x_coho1_skagit_wo_hatchery$pred$.rownames), 8,12)

skagit_covariates_coho1 <- as.data.frame(t(covariates_coho1_skagit_night))

glimpse(skagit_covariates_coho1)

skagit_covariates_coho1_long <- skagit_covariates_coho1 %>% 
  mutate(doy = as.numeric(rownames(skagit_covariates_coho1))) %>%
  pivot_longer(cols = -c(doy), 
               names_to = c(".value","year", "daynight_category","trap"),
               names_pattern = "(.*)_(.{4})_(.{5})_(.{5})")
head(skagit_covariates_coho1_long)

#merge skagit_covariates_coho1_long$coho1_hatchery_perhour_inp with predict_coho1_skagit$pred

predict_x_coho1_skagit_wo_hatchery$pred2 <- predict_x_coho1_skagit_wo_hatchery$pred %>% 
  mutate(doy = t+100) %>% 
  left_join(skagit_covariates_coho1_long, by = c("year", "daynight_category", "trap","doy"))

skagit_covariates_coho1$doy <- as.numeric(rownames(skagit_covariates_coho1))


ggplot(data = predict_x_coho1_skagit_wo_hatchery$pred %>% 
         filter(year == 2021 & daynight_category == "night" & trap == 'scoop'))+
  
  geom_line(aes(x = t+100, y = estimate, color = "predicted"))+
  geom_point(aes(x = t+100, y = .x, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = skagit_covariates_coho1, aes(x = doy, y = percent_hatchery_diff_2016_night_scoop), 
            color = "cadetblue", linewidth = 1.2, alpha = 0.7) +
  scale_color_manual(name = "Legend", values = c("predicted" = "salmon", "observed" = "salmon"))+
  labs(x = "Time", y = "log(scaled(coho salmon per hour) + 1)", title = "")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.2, linetype="solid"
    ),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )


#figure out the scale value


data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In


data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                       na.rm = FALSE)

lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)


data_day_night <- data_day_night %>%
  filter(doy > 100, doy < 200) %>%
  group_by(year, daytime_category) %>% 
  mutate(coho1_wild_perhour_cumsum = cumsum(ifelse(is.na(coho1_wild_perhour),
                                                   0,
                                                   coho1_wild_perhour))) %>% 
  mutate(coho1_proportion = coho1_wild_perhour_cumsum/
           sum(coho1_wild_perhour,na.rm = TRUE))

median <- data_day_night %>% 
  group_by(year, daytime_category) %>% 
  dplyr::select(year,coho1_proportion,doy) %>%
  summarise(median_doy = doy[which.min(abs(coho1_proportion - 0.5))])

data_day_night <- left_join(data_day_night,median,by = c("year","daytime_category"))



#calculate the difference between the median day of migration and the day of year

data_day_night <- data_day_night %>%
  mutate(season = median_doy - doy)


subset_coho_summer_perhour_night <- arrange(data_day_night,doy) %>%
  filter(doy > 100 & doy < 150 & daytime_category == 'night') %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

# make data frame to store mean and sd of each row and year and trap

mean_sd <- data.frame(matrix(ncol = 3, nrow = dim(subset_coho_summer_perhour_night)[1]))


for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
  # print(sd(subset_coho_summer_perhour_night[i,], na.rm = TRUE))
  # print(scale(subset_coho_summer_perhour_night[i,]))
  # subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
  mean_sd[i,1] = mean(subset_coho_summer_perhour_night[i,], na.rm = TRUE)
  mean_sd[i,2] = sd(subset_coho_summer_perhour_night[i,], na.rm = TRUE)
  mean_sd[i,3] = substr(rownames(subset_coho_summer_perhour_night)[i], 1, 4)
  mean_sd[i,4] = substr(rownames(subset_coho_summer_perhour_night)[i], 12, 17)
}

colnames(mean_sd) <- c("mean", "sd", "year", "trap")

#take the predicted values, multiply by the scale value, and then add the mean value,
#then exponentiate to get the original value

predict_x_coho1_skagit_wo_hatchery$pred3 <- predict_x_coho1_skagit_wo_hatchery$pred2 %>% 
  left_join(mean_sd, by = c("year", "trap")) %>%
  mutate(estimate_original = exp((estimate * sd) + mean) - 1,
         .x_original = exp((.x * sd) + mean) - 1,
         `Lo 95_original` = exp((`Lo 95` * sd) + mean) - 1,
         `Hi 95_original` = exp((`Hi 95` * sd) + mean) - 1)

predict_x_coho1_skagit_w_hatchery$pred3 <- predict_x_coho1_skagit_w_hatchery$pred2 %>% 
  left_join(mean_sd, by = c("year", "trap")) %>%
  mutate(estimate_original = exp((estimate * sd) + mean) - 1,
         .x_original = exp((.x * sd) + mean) - 1,
         `Lo 95_original` = exp((`Lo 95` * sd) + mean) - 1,
         `Hi 95_original` = exp((`Hi 95` * sd) + mean) - 1)

w_hatchery <- as.data.frame(predict_x_coho1_skagit_w_hatchery$pred3) %>% 
  filter(year == 2016 & daynight_category == "night" & trap == 'scoop',doy == 125) %>%
  dplyr::select(estimate_original)

wo_hatchery <- as.data.frame(predict_x_coho1_skagit_wo_hatchery$pred3) %>%
  filter(year == 2016 & daynight_category == "night" & trap == 'scoop',doy == 125) %>%
  dplyr::select(estimate_original)

p1 <- ggplot(data = predict_x_coho1_skagit_wo_hatchery$pred3 %>% 
               filter(year == 2016 & daynight_category == "night" & trap == 'scoop'))+
  
  geom_line(aes(x = t+100, y = estimate_original, color = "predicted wo hatchery"), linetype = "dashed")+
  geom_point(aes(x = t+100, y = .x_original, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95_original`, ymax = `Hi 95_original`), alpha = 0.1)+
  # geom_line(data = skagit_covariates_coho1, aes(x = doy, y = percent_hatchery_diff_2016_night_scoop),#*4.493336),
  #           color = "cadetblue", linewidth = 1.2, alpha = 0.7) +
  geom_line(data = data_day_night %>%
              filter(year == 2016 & daytime_category == "night" & trap == 'scoop',
                     doy > 100 & doy < 150),
            aes(x = doy, y = coho1_hatchery_perhour_inp,color = "hatchery"),
            linewidth = 1.2, alpha = 0.7
  )+
  geom_line(data = predict_x_coho1_skagit_w_hatchery$pred3 %>% 
              filter(year == 2016 & daynight_category == "night" & trap == 'scoop'),
            aes(x = t+100, y = estimate_original, color = "predicted"))+
  # geom_point(data = predict_x_coho1_skagit_w_hatchery$pred %>% 
  #             filter(year == 2016 & daynight_category == "night" & trap == 'scoop'),
  #           aes(x = t+100, y = .x_original, color = "observed"), alpha = 0.5)+
  geom_ribbon(data = predict_x_coho1_skagit_w_hatchery$pred3 %>% 
                filter(year == 2016 & daynight_category == "night" & trap == 'scoop'),
              aes(x = t+100, ymin = `Lo 95_original`, ymax = `Hi 95_original`), alpha = 0.1)+
  
  
  scale_color_manual(name = "Legend", values = c("predicted" = "salmon", 
                                                 "observed" = "salmon",
                                                 "predicted wo hatchery" = "slategray",
                                                 "hatchery" = "cadetblue"
                                                 ))+
  
  
  geom_text(data = w_hatchery, aes(x = 125, y = estimate_original, label = round(estimate_original, 2)), 
            hjust = +1.1, vjust = 0, size = 3, color = "salmon")+
  geom_text(data = wo_hatchery, aes(x = 125, y = estimate_original, label = round(estimate_original, 2)), 
            hjust = +1.2, vjust = 0, size = 3, color = "slategray")+
  
  labs(x = "Time", y = "Number of coho salmon per hour", title = "2016, scoop")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.2, linetype="solid"
    ),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  )


p1

w_hatchery <- as.data.frame(predict_x_coho1_skagit_w_hatchery$pred3) %>% 
  filter(year == 2020 & daynight_category == "night" & trap == 'scoop',doy == 125) %>%
  dplyr::select(estimate_original)

wo_hatchery <- as.data.frame(predict_x_coho1_skagit_wo_hatchery$pred3) %>%
  filter(year == 2020 & daynight_category == "night" & trap == 'scoop',doy == 125) %>%
  dplyr::select(estimate_original)

p2 <- ggplot(data = predict_x_coho1_skagit_wo_hatchery$pred3 %>% 
               filter(year == 2020 & daynight_category == "night" & trap == 'scoop'))+
  
  geom_line(aes(x = t+100, y = estimate_original, color = "predicted \nwithout hatchery" ), linetype = "dashed")+
  geom_point(aes(x = t+100, y = .x_original, color = "observed"), alpha = 0.5)+
  geom_ribbon(aes(x = t+100, ymin = `Lo 95_original`, ymax = `Hi 95_original`), alpha = 0.2)+
  
  geom_line(data = data_day_night %>%
              filter(year == 2020 & daytime_category == "night" & trap == 'scoop',
                     doy > 100 & doy < 150),
            aes(x = doy, y = coho1_hatchery_perhour_inp,color = "hatchery"),
            linewidth = 1.2, alpha = 0.7
  )+
  # geom_line(data = skagit_covariates_coho1, aes(x = doy, y = percent_hatchery_diff_2020_night_scoop,color = "hatchery"),#*4.493336),
  #           linewidth = 1.2, alpha = 0.7)+
  geom_line(data = predict_x_coho1_skagit_w_hatchery$pred3 %>% 
              filter(year == 2020 & daynight_category == "night" & trap == 'scoop'),
            aes(x = t+100, y = estimate_original, color = "predicted \nwith hatchery"))+
  
  geom_ribbon(data = predict_x_coho1_skagit_w_hatchery$pred3 %>% 
                filter(year == 2020 & daynight_category == "night" & trap == 'scoop'),
              aes(x = t+100, ymin = `Lo 95_original`, ymax = `Hi 95_original`), alpha = 0.2)+
  
  
  scale_color_manual(name = "", values = c("predicted \nwith hatchery" = "salmon", 
                                           "observed" = "salmon",
                                           "predicted \nwithout hatchery" = "slategray",
                                           "hatchery" = "cadetblue"
  ))+
  
  geom_text(data = w_hatchery, aes(x = 120, y = estimate_original, label = round(estimate_original,2)), 
            color = "salmon", size = 3, vjust = 0)+
  geom_text(data = wo_hatchery, 
            aes(x = 120, y = estimate_original, label = round(estimate_original,2)),
            color = "slategray", size = 3, vjust = 0)+
  
  
  labs(x = "Time", y = "Number of coho salmon per hour", title = "2020, scoop")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=1.2, linetype="solid"
    ),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "right"
  )

p2
library(patchwork)
p<- p1+p2
ggsave(here("visualizations","output","predictions_w_wo_hatchery.png"), width = 10, height = 5)


#caluclate sd of each of the rows of the covariates 
# save in dataframe with year, trap, sd, covariate_name

covariates_coho1_skagit_night_temp <- arrange(data_day_night,doy) %>%
  filter(doy >100 & doy < 150 & daytime_category == 'night') %>%
  dplyr::select(year,doy, daytime_category, temp, flow, photoperiod, atu_solstice,
                lunar_phase, resid, temp_diff, flow_diff, photo_diff,
                coho1_hatchery_perhour_inp, trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = c(
    temp, flow, photoperiod, atu_solstice, lunar_phase,
    resid, temp_diff, flow_diff, photo_diff, coho1_hatchery_perhour_inp), names_vary = "fastest") %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

sd_df <- data.frame(matrix(NA, nrow = 24, 
                           ncol = 4))
colnames(sd_df) <- c("year", "trap", "covariate", "sd")

for(i in 217:nrow(covariates_coho1_skagit_night_temp)){
  #year is between the first two '_' of the rownames
  sd_df[i-216,1] <- strsplit(rownames(covariates_coho1_skagit_night_temp)[i], "_")[[1]][5]
  sd_df[i-216,2] <- strsplit(rownames(covariates_coho1_skagit_night_temp)[i], "_")[[1]][6]
  sd_df[i-216,3] <- strsplit(rownames(covariates_coho1_skagit_night_temp)[i], "_")[[1]][7]
  sd_df[i-216,4] <- sd(covariates_coho1_skagit_night_temp[i,])
}

#add this to the covariates_coho1_skagit_night_temp dataframe





(w_hatchery - wo_hatchery)*100/wo_hatchery

