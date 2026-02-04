# try best model wo 2019 data


# Load libraries

library(MARSS)
library(broom)
library(here)
library(zoo)
library(MASS)
library(modelr)
library(qpcR)
library(GGally)
library(ggplot2)
library(tidyverse)

# Load data

#load data
data <- read.csv(here("puyallup", "data","puyallup_final.csv"),
                 na.strings = c("NA",""))

data <- data %>%
  filter(doy > 90, doy < 230) %>%
  group_by(year) %>% 
  mutate(chinook0_wild_perhour_day_cumsum = cumsum(ifelse(is.na(chinook0_wild_perhour_day),
                                                          0,
                                                          chinook0_wild_perhour_day)),
         chinook0_wild_perhour_night_cumsum = cumsum(ifelse(is.na(chinook0_wild_perhour_night),
                                                            0,
                                                            chinook0_wild_perhour_night))) %>% 
  mutate(chinook0_day_proportion = chinook0_wild_perhour_day_cumsum/
           sum(chinook0_wild_perhour_day,na.rm = TRUE),
         chinook0_night_proportion = chinook0_wild_perhour_night_cumsum/
           sum(chinook0_wild_perhour_night,na.rm = TRUE))

data <- data %>% 
  group_by(year) %>% 
  mutate(chinook0_hatchery_perhour_day_diff = c(NA,diff(chinook0_hatchery_perhour_day,1))) %>% 
  mutate(chinook0_hatchery_perhour_night_diff = c(NA,diff(chinook0_hatchery_perhour_night,1))) %>%
  mutate(coho1_hatchery_perhour_day_diff = c(NA,diff(coho1_hatchery_perhour_day,1))) %>%
  mutate(coho1_hatchery_perhour_night_diff = c(NA,diff(coho1_hatchery_perhour_night,1))) %>%
  mutate(temp_day_rolling_mean = rollmean(temp_day, k = 32, fill = NA, align = "right"),
         temp_night_rolling_mean = rollmean(temp_night, k = 32, fill = NA, align = "right"),
         flow_day_rolling_mean = rollmean(flow_day, k = 32, fill = NA, align = "right"),
         flow_night_rolling_mean = rollmean(flow_night, k = 32, fill = NA, align = "right"),
         temp_anomaly_day = temp_day - temp_day_rolling_mean,
         temp_anomaly_night = temp_night - temp_night_rolling_mean,
         flow_anomaly_day = flow_day - flow_day_rolling_mean,
         flow_anomaly_night = flow_night - flow_night_rolling_mean)

#calculate the median day of migration for each year (day of year that 0.5 proportion of fish have passed)

# make column for median day of migration - day of year and call it season

median <- data %>% 
  group_by(year) %>% 
  dplyr::select(year,chinook0_day_proportion,
                chinook0_night_proportion,
                doy) %>%
  summarise(median_day_doy = doy[which.min(abs(chinook0_day_proportion - 0.5))],
            median_night_doy = doy[which.min(abs(chinook0_night_proportion - 0.5))])


data <- left_join(data,median,by = c("year"))

#calculate the difference between the median day of migration and the day of year

data <- data %>%
  mutate(season_day = median_day_doy - doy,
         season_night = median_night_doy - doy)


covariates_chinook0_puyallup_w_temp <- arrange(data,doy) %>%
  filter(year!=2019) %>% 
  filter(doy >130 & doy <= 218) %>%
  dplyr::select(year, doy, flow_anomaly_day, flow_anomaly_night, 
                temp_anomaly_day, temp_anomaly_night,
                # secchi_depth_day, secchi_depth_night,
                # lunar_phase_day, lunar_phase_night, 
                season_day, season_night,
                flow_diff_day, flow_diff_night, 
                # photo_diff_day, photo_diff_night, 
                temp_diff_day, temp_diff_night, 
                # resid_day, resid_night,
                chinook0_hatchery_perhour_day_diff, chinook0_hatchery_perhour_night_diff) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_anomaly_day, flow_anomaly_night, 
    temp_anomaly_day, temp_anomaly_night,
    # secchi_depth_day, secchi_depth_night,
    # lunar_phase_day, lunar_phase_night,
    season_day, season_night, 
    flow_diff_day, flow_diff_night, 
    # photo_diff_day, photo_diff_night, 
    temp_diff_day, temp_diff_night,
    # resid_day, resid_night,
    chinook0_hatchery_perhour_day_diff, chinook0_hatchery_perhour_night_diff)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()




#scaling the variables

num_years = 2021-2004+1-1
num_rows = num_years*2
total_covariates = dim(covariates_chinook0_puyallup_w_temp)[1]


for(i in 1:(num_rows*2)){ # everything except diffs and hatchery
  # print(rownames(covariates_chinook0_puyallup_w_temp)[i])
  covariates_chinook0_puyallup_w_temp[i,] = scale(covariates_chinook0_puyallup_w_temp[i,])[,1]
}

#just scale

for(i in (num_rows*2 + 1):(total_covariates)){
  # print(rownames(covariates_chinook0_puyallup_w_temp)[i])
  covariates_chinook0_puyallup_w_temp[i,] = scale(covariates_chinook0_puyallup_w_temp[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_chinook_summer_perhour <- arrange(data,doy) %>%
  filter(year!=2019) %>% 
  filter(doy > 130 & doy <= 218) %>%
  mutate(log.value_day = log(chinook0_wild_perhour_day + 1), 
         log.value_night = log(chinook0_wild_perhour_night + 1)) %>%
  dplyr::select(log.value_day, log.value_night ,year,doy) %>%
  pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour)[1]){
  subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
}


num_rows = num_years*2

c<-NULL

for(kk in c(1,2,3,4,6)){
  c = rbind(c,covariates_chinook0_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),])
  name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-5)
  print(name_individual)
  
  # out=data.frame(c=name_individual, d = "None",
  #                logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
  #                num.iter=fit$numIter, converged=!fit$convergence,
  #                stringsAsFactors = FALSE)
  # out.tab.all.years_night=rbind(out.tab.all.years_night,out)
  # fits.all.years_night=c(fits.all.years_night,list(fit))
}
c_num = 5 # num covariates
fit.model = c(list(c= c), mod_list(num_rows,c_num,1, FALSE, FALSE))

fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

ci <- tidy(fit)
puyallup_chinook_plot_anomaly_wo_2019 <- ggplot(ci[c(37:42),], 
                                        aes(x = c("Flow\n anomaly",
                                                  "Temperature\n anomaly",
                                                  "Season", 
                                                  "Flow\n difference",
                                                  "Hatchery\ndifference, day", 
                                                  "Hatchery\ndifference, night"),
                                            y = estimate, 
                                            ymin = conf.low, 
                                            ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Puyallup River, Chinook sub yearlings"
  )+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits = c("Flow\n anomaly","Season","Temperature\n anomaly","Flow\n difference",
               "Hatchery\ndifference, day","Hatchery\ndifference, night"
               
    )) + 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  coord_flip()+
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.35,0.35))+
  geom_rect(aes(xmin = 4.5, xmax = 6.5, ymin = -0.34, ymax = 0.34), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 0.7, y = 0.3, label = "Puyallup"), size = 10)

puyallup_chinook_plot_anomaly_wo_2019

ggsave(here("visualizations",
            "output",
            "puyallup_chinook_covariates_estimates_anomaly_wo_2019.jpeg"), width = 16, height = 6)

