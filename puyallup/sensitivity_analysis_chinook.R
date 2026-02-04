
# load libraries

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
#sensitivity analysis

#best model had all covariates except temp difference

#loop through different percent values
percent <- c(0,0.05,0.1,0.15)

sensitivity <- data.frame(matrix(NA, nrow = 6*length(percent), ncol = 3))

sensitivity$variable <- rep(c("flow anomaly","temperature anomaly",
                              "season", "flow difference",
                              "hatchery difference day",
                              "hatchery difference night"),length(percent))

sensitivity$percent <- rep(percent, each = 6)


for(p in 1:length(percent)){
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
           flow_anomaly_night = flow_night - flow_night_rolling_mean) %>% 
    
    mutate(percent_wild_day = chinook0_wild_perhour_day + percent[p]*chinook0_hatchery_perhour_day,
           percent_wild_night = chinook0_wild_perhour_night + percent[p]*chinook0_hatchery_perhour_night,
           percent_hatchery_day = chinook0_hatchery_perhour_day - percent[p]*chinook0_hatchery_perhour_day,
           percent_hatchery_night = chinook0_hatchery_perhour_night - percent[p]*chinook0_hatchery_perhour_night) %>%
    mutate(percent_hatchery_day_diff = c(NA, diff(percent_hatchery_day,1)),
           percent_hatchery_night_diff = c(NA, diff(percent_hatchery_night,1)),
           chinook0_percent_wild_perhour_day_cumsum = cumsum(ifelse(is.na(percent_wild_day),
                                                             0,
                                                             percent_wild_day)),
           chinook0_percent_wild_perhour_night_cumsum = cumsum(ifelse(is.na(percent_wild_night),
                                                                    0,
                                                                    percent_wild_night))) %>%
    mutate(chinook0_percent_wild_day_proportion = chinook0_percent_wild_perhour_day_cumsum/
             sum(percent_wild_day,na.rm = TRUE),
           chinook0_percent_wild_night_proportion = chinook0_percent_wild_perhour_night_cumsum/
             sum(percent_wild_night,na.rm = TRUE))
  
  
  median <- data %>% 
    group_by(year) %>% 
    dplyr::select(year,chinook0_percent_wild_day_proportion,
                  chinook0_percent_wild_night_proportion,
                  doy) %>%
    summarise(median_percent_day_doy = doy[which.min(abs(chinook0_percent_wild_day_proportion - 0.5))],
              median_percent_night_doy = doy[which.min(abs(chinook0_percent_wild_night_proportion - 0.5))])
  
  
  data <- left_join(data,median,by = c("year"))
  
  #calculate the difference between the median day of migration and the day of year
  
  data <- data %>%
    mutate(season_percent_day = median_percent_day_doy - doy,
           season_percent_night = median_percent_night_doy - doy)
  
  
  covariates_chinook0_puyallup_w_temp <- arrange(data,doy) %>%
    filter(doy >130 & doy <= 218) %>%
    dplyr::select(year, doy, flow_anomaly_day, flow_anomaly_night, 
                  temp_anomaly_day, temp_anomaly_night,
                  # secchi_depth_day, secchi_depth_night,
                  # lunar_phase_day, lunar_phase_night, 
                  season_percent_day, season_percent_night,
                  flow_diff_day, flow_diff_night, 
                  # photo_diff_day, photo_diff_night, 
                  temp_diff_day, temp_diff_night, 
                  # resid_day, resid_night,
                  percent_hatchery_day_diff, percent_hatchery_night_diff) %>%
    pivot_wider(names_from = c(year), values_from = c(
      flow_anomaly_day, flow_anomaly_night, 
      temp_anomaly_day, temp_anomaly_night,
      # secchi_depth_day, secchi_depth_night,
      # lunar_phase_day, lunar_phase_night,
      season_percent_day, season_percent_night, 
      flow_diff_day, flow_diff_night, 
      # photo_diff_day, photo_diff_night, 
      temp_diff_day, temp_diff_night,
      # resid_day, resid_night,
      percent_hatchery_day_diff, percent_hatchery_night_diff)) %>%
    column_to_rownames(var = "doy") %>%
    as.matrix() %>%
    t()
  
  
  #scaling the variables
  
  num_years = 2021-2004+1
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
    filter(doy > 130 & doy <= 218) %>%
    mutate(log.value_day = log(percent_wild_day + 1), 
           log.value_night = log(percent_wild_night + 1)) %>%
    dplyr::select(log.value_day, log.value_night ,year,doy) %>%
    pivot_wider(names_from = c(year), values_from = c(log.value_day, log.value_night)) %>%
    column_to_rownames(var = "doy") %>%
    as.matrix() %>%
    t()
  
  for(i in 1:dim(subset_chinook_summer_perhour)[1]){
    subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
  }
  
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
  
  sensitivity[(1+6*(p-1)):(6*p),1] <- ci$estimate[39:44]
  sensitivity[(1+6*(p-1)):(6*p),2] <- ci$conf.low[39:44]
  sensitivity[(1+6*(p-1)):(6*p),3] <- ci$conf.up[39:44]
  
}

#rename X1 as estimate
colnames(sensitivity)[1:4] <- c("estimate", "conf.low", "conf.up", "variable")


p5 <- sensitivity %>% 
  filter(variable == "hatchery difference day" | variable == "hatchery difference night") %>%
  ggplot(aes(x = percent, y = estimate, group = variable, color = variable)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up)) +
  geom_line(alpha = 0.8) +
  xlim(-0.2,0.2) +
  # ylim(-0.1,0.1) +
  stat_smooth(method="lm",fullrange=TRUE, alpha = 0.5, aes(fill = variable)) +
  # geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#6e90a5", "gray30"),
                     labels = c("Hatchery difference day", "Hatchery difference night")) +
  scale_fill_manual(values = c("#6e90a5", "gray30"),
                    labels = c("Hatchery difference day", "Hatchery difference night")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(#title = "Sensitivity Analysis - Skagit River Coho Salmon",
    x = "Proportion of hatchery chinook relabelled as wild chinook",
    y = "Estimate of effect of hatchery difference") +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank())

p5

ggsave(here("puyallup", "output", "sensitivity_analysis_chinook.png"),
       plot = p5,
       width = 6,
       height = 4,
       dpi = 300)



