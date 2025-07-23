# Goal - To run marss analysis on chinook data from Puyallup River with
# median day of migration - day of year for every year as covariate
# to do model selection with other covariates

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
  filter(doy > 100, doy < 230) %>%
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



#calculate the median day of migration for each year (day of year that 0.5 proportion of fish have passed)

# make column for median day of migration - day of year and call it season

median <- data %>% 
  group_by(year) %>% 
  select(year,chinook0_day_proportion,
         chinook0_night_proportion,
         doy) %>%
  summarise(median_day_doy = doy[which.min(abs(chinook0_day_proportion - 0.5))],
            median_night_doy = doy[which.min(abs(chinook0_night_proportion - 0.5))])
            

data <- left_join(data,median,by = c("year"))

#calculate the difference between the median day of migration and the day of year

data <- data %>%
  mutate(season_day = median_day_doy - doy,
         season_night = median_night_doy - doy)

data %>% 
  ungroup() %>%
  filter(doy > 130 & doy <= 218) %>%
  select(flow_day,lunar_phase_day, season_day, resid_day, temp_diff_day, flow_diff_day, 
         photo_diff_day,secchi_depth_day) %>%
  GGally::ggpairs(aes(alpha = 0.2))

ggsave(here("puyallup","output","chinook_covariates_correlation.png"),width = 10, height = 10)


data %>%
  filter(doy > 130 & doy <= 218, year == 2019) %>%
  ggplot() +
  # geom_line(aes(x = doy, y = scale(photoperiod), group = year), col = "darkred") +
  geom_line(aes(x = doy, y = scale(lunar_phase_night), group = year), col = "darkblue") +
  geom_line(aes(x = doy, y = log(chinook0_wild_perhour_night+1), group = year), col = "darkgreen") +
  # geom_line(aes(x = doy, y = scale(season, center = FALSE), group = year), col = "darkred") +
  # ylim(-1,1)+
  # geom_line(aes(x = doy, y = photo_diff, group = year), col = "darkblue") +
  # geom_line(aes(x = doy, y = temp_diff, group = year), col = "darkgreen") +
  # geom_line(aes(x = doy, y = flow_diff, group = year), col = "darkorange") +
  # geom_line(aes(x = doy, y = flow, group = year), col = "darkviolet") +
  # geom_line(aes(x = doy, y = lunar_phase, group = year), col = "darkgrey") +
  # geom_line(aes(x = doy, y = resid, group = year), col = "black") +
  theme_minimal()


#plot temperature residuals for each year as a function of day of year (doy) using ggplot2

data %>% 
  filter(doy>90, doy <= 218) %>% 
  ggplot() +
  geom_line(aes(x = doy, y = resid_day, group = as.factor(year), col = "residuals"),alpha = 0.5) +
  geom_line(aes(x = doy, y = temp_day, group= as.factor(year), col = "temperature"),alpha = 0.5)+
  geom_line(aes(x = doy, y = photoperiod_day, group = as.factor(year), col = "photoperiod"), 
            linewidth = 2,alpha = 0.7)+
  scale_color_manual(name = "",values = c("residuals" = "slategray",
                                "temperature" = "#B17A79",
                                "photoperiod" = "#E8D68A")) +
  theme_classic() +
  labs(x = "Day of Year",
       y = "Value",
       title = "") +
  theme(legend.position = "right")

ggsave(here("puyallup","output","puyallup_temperature_residuals.png"),width = 6, height = 5)


covariates_chinook0_puyallup_w_temp <- arrange(data,doy) %>%
  filter(doy >130 & doy <= 218) %>%
  dplyr::select(year,doy, flow_day, flow_night, 
                secchi_depth_day, secchi_depth_night,
                lunar_phase_day, lunar_phase_night, 
                season_day, season_night,
                flow_diff_day, flow_diff_night, 
                photo_diff_day, photo_diff_night, 
                temp_diff_day, temp_diff_night, 
                resid_day, resid_night,
                chinook0_hatchery_perhour_day, chinook0_hatchery_perhour_night) %>%
  pivot_wider(names_from = c(year), values_from = c(
    flow_day, flow_night, 
    secchi_depth_day, secchi_depth_night,
    lunar_phase_day, lunar_phase_night,
    season_day, season_night, 
    flow_diff_day, flow_diff_night, photo_diff_day, photo_diff_night, 
    temp_diff_day, temp_diff_night,
    resid_day, resid_night,
    chinook0_hatchery_perhour_day, chinook0_hatchery_perhour_night)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()




#scaling the variables

num_years = 2021-2004+1
num_rows = num_years*2
total_covariates = dim(covariates_chinook0_puyallup_w_temp)[1]


for(i in 1:(num_rows*3)){ # everything except diffs and hatchery
  # print(rownames(covariates_chinook0_puyallup_w_temp)[i])
  covariates_chinook0_puyallup_w_temp[i,] = scale(covariates_chinook0_puyallup_w_temp[i,])[,1]
}

#just scale

for(i in (num_rows*3 + 1):(total_covariates)){
  # print(rownames(covariates_chinook0_puyallup_w_temp)[i])
  covariates_chinook0_puyallup_w_temp[i,] = scale(covariates_chinook0_puyallup_w_temp[i,], center = FALSE, scale= TRUE)[,1]
}



#subset response variable
subset_chinook_summer_perhour <- arrange(data,doy) %>%
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


#function for C matrix
#######


Cmat <- function(nyears,ncov,hatchery=0, day_on_night = FALSE){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  
  if(hatchery == 1){
    C <- matrix(list(0),nrow = nyears,ncol = ncov*nyears, byrow = TRUE)
    for(i in 1:nyears){
      for(j in 1:(ncov*nyears)){
        
        for(k in 1:ncov){
          # print(((k-1)*nyears+i))
          if(j == ((k-1)*nyears+i)){
            # print(i)
            # print(j)
            # print(vars[k])
            if(k==ncov){
              if(i<=nyears/2){
                C[i,j] <- "day"
                if(day_on_night){
                  C[i+nyears/2,j] <- "day_on_night"
                }
                
              }
              else{
                C[i,j] <- "night"
              }
            }
            
            else{
              C[i,j] <- vars[k]
            }
            
            
          }
          
        }
      }
      
    }
    
  }
  else{
    C <- matrix(list(0),nrow = nyears,ncol = ncov*nyears, byrow = TRUE)
    for(i in 1:nyears){
      for(j in 1:(ncov*nyears)){
        
        for(k in 1:ncov){
          # print(((k-1)*nyears+i))
          if(j == ((k-1)*nyears+i)){
            # print(i)
            # print(j)
            # print(vars[k])
            C[i,j] <- vars[k]
            
          }
          
        }
      }
      
    }
  }
  
  return(C)
}


#######

#function for Q matrix
#######


Qmat <- function(nyears){
  Q <- matrix(list(0),nrow = nyears,ncol = nyears, byrow = TRUE)
  for(i in 1:nyears){
    for(j in 1:nyears){
      if(i==j){
        if(i <= nyears/2){
          Q[i,j] <- "q_d"
        }
        else{
          Q[i,j] <- "q_n"
        }
      }
    }
  }
  return(Q)
}


#######

#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0, day_on_night = FALSE, unequal_q = FALSE){
  
  if(unequal_q){
    Q = Qmat(nyears)
    
  }
  else{
    Q = "diagonal and equal"
  }
  
  if(ncov == 0){
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = Q
    )
  }
  else{
    if((ncov == 1) & (hatchery == 0)){
      C = "diagonal and equal"
    }
    else{
      C = Cmat(nyears,ncov,hatchery, day_on_night)
    }
    mod.list = list(
      B = "identity",
      U = "zero",
      Z = "identity",
      A = "zero",
      R = "diagonal and equal",
      Q = Q,
      C = C
    )
  }
  
  return(mod.list)
}

######

#function to get covariate combinations
##### 
get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}

nyears = num_years*2
c = NULL
fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,FALSE))
fit_equal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                     control=list(maxit=2000))
fit_equal_q$AICc

fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,TRUE))
fit_unequal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                       control=list(maxit=2000))
fit_unequal_q$AICc

df_errors <- data.frame(Error_structure = c("Equal","Unequal"), 
                        AICc = c(fit_equal_q$AICc,fit_unequal_q$AICc))

#save table
write.csv(df_errors, file = here("puyallup",
                                 "output","error_structure_puyallup_chinook.csv"))


out.tab.all.years <- NULL
fits.all.years <- NULL
nyears = num_years*2
for(kk in c(1,2,4,6,8)){
  c = covariates_chinook0_puyallup_w_temp[((1+(kk-1)*num_rows):(kk*num_rows)),]
  name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(kk-1)*num_rows]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0, FALSE, FALSE))
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  ci = tidy(fit)
  out=data.frame(c=name_individual, Estimate = ci[39,2],
                 conf_low = ci[39,4], conf_high = ci[39,5],
                 AICc=fit$AICc)
  out.tab.all.years=rbind(out.tab.all.years,out)
  fits.all.years=c(fits.all.years,list(fit))
}

out.tab.all.years$delta_AICc = out.tab.all.years$AICc - min(out.tab.all.years$AICc)
out.tab.all.years_ordered = out.tab.all.years[order(out.tab.all.years$delta_AICc),]
out.tab.all.years_ordered

#save table
write.csv(out.tab.all.years_ordered, file = here("puyallup","output",
                                                 "chinook_correlated_covariates_selection.csv"))
#will go with season, it has delta aicc on 0.17

#flow, lunar phase, season, flow_diff, temp_diff, hatchery

num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:6)
out.tab_season_puyallup <- NULL
out_riv_puyallup <- NULL
fits_season_puyallup <- list()

amount_done = length(fits_season_puyallup)
for(i in (amount_done + 1):length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  season = 0
  temperature_difference = 0
  flow = 0
  lunar_phase = 0
  hatchery = 0
  flow_difference = 0
  for(j in covariates){
    if(j == 1){
      k = 1
      flow =1
    }
    else if(j==2){
      k = 3
      lunar_phase = 1
    }
    else if(j==3){
      k = 4
      season = 1
    }
    else if(j==4){
      k = 7
      temperature_difference = 1
    }
    
    else if(j==5){
      k = 5
      flow_difference = 1
    }
    else if(j==6){
      k = 9
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_puyallup_w_temp[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==9){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, FALSE))
  }else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, FALSE))
  }
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season,
                 temperature_difference = temperature_difference, flow = flow,
                 flow_difference = flow_difference,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv_puyallup <- rbind(out,out_riv_puyallup)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    temperature_difference = ifelse(temperature_difference==1, fit$par$U[which(covariates == 4),], NA),
                    flow = ifelse(flow == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_difference = ifelse(flow_difference == 1, fit$par$U[which(covariates == 5),], NA),
                    lunar_phase = ifelse(lunar_phase == 1, fit$par$U[which(covariates == 2),], NA),
                    hatchery_day = ifelse(hatchery == 1, fit$par$U["day",], NA),
                    hatchery_night = ifelse(hatchery == 1, fit$par$U["night",], NA),
                    AICc=fit$AICc)
  out.tab_season_puyallup=rbind(out.tab_season_puyallup,out2)
  fits_season_puyallup=c(fits_season_puyallup,list(fit))
  
  
}

out.tab_season_puyallup$deltaAICc <- out.tab_season_puyallup$AICc - min(out.tab_season_puyallup$AICc)
min.AICc <- order(out.tab_season_puyallup$AICc)
out.tab_season_puyallup.ordered <- out.tab_season_puyallup[min.AICc, ]
out.tab_season_puyallup.ordered

#round all the values to 2 decimal places

out.tab_season_puyallup.ordered <- round(out.tab_season_puyallup.ordered, 2)

write.csv(out.tab_season_puyallup.ordered, 
          file = here("puyallup",
                      "output","model_selection_puyallup_chinook.csv"))

#best model

puyallup_chinook_best_model <- fits_season_puyallup[[which.min(sapply(fits_season_puyallup, function(x) x$AICc))]]

#save

save(puyallup_chinook_best_model, file = here("puyallup",
                                               "output","puyallup_chinook_best_model.RData"))

#ci
tidy(puyallup_chinook_best_model)

autoplot(puyallup_chinook_best_model)


#relative importance

fit.model = c(mod_list(num_rows,0,0, FALSE, FALSE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0


out=data.frame(c=name, season = season,
               temperature_difference = temperature_difference, flow = flow,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out_riv_puyallup$deltaAICc <- NULL
out_riv_puyallup <- rbind(out,out_riv_puyallup)


out_riv_puyallup$rel.LL <- NULL
out_riv_puyallup$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_season_puyallup=c(fits_season_puyallup,list(fit))

weights <- akaike.weights(out_riv_puyallup$AICc)

out_riv_puyallup$deltaAICc <- weights$deltaAIC
out_riv_puyallup$rel.LL <- weights$rel.LL
out_riv_puyallup$weights <- weights$weights


min.AICc <- order(out_riv_puyallup$AICc)
out_riv_puyallup.ordered <- out_riv_puyallup[min.AICc, ]
out_riv_puyallup.ordered

out_riv_puyallup.ordered$cumulative_weights <- cumsum(out_riv_puyallup.ordered$weights)

relative_importance_season <- sum(out_riv_puyallup$weights[out_riv_puyallup$season==1])
relative_importance_temperature_difference <- sum(out_riv_puyallup$weights[out_riv_puyallup$temperature_difference==1])
relative_importance_flow <- sum(out_riv_puyallup$weights[out_riv_puyallup$flow==1])
relative_importance_flow_difference <- sum(out_riv_puyallup$weights[out_riv_puyallup$flow_difference==1])
relative_importance_lunar_phase <- sum(out_riv_puyallup$weights[out_riv_puyallup$lunar_phase==1])
relative_importance_hatchery <- sum(out_riv_puyallup$weights[out_riv_puyallup$hatchery==1])

riv_season_puyallup <- data.frame(variable = c("season",
                                      "temperature difference",
                                      "flow",
                                      "flow difference",
                                      "lunar phase",
                                      "hatchery"),
                         relative_importance = c(relative_importance_season,
                                                 relative_importance_temperature_difference,
                                                 relative_importance_flow,
                                                 relative_importance_flow_difference,
                                                 relative_importance_lunar_phase,
                                                 relative_importance_hatchery))


#save
write.csv(riv_season_puyallup[order(riv_season_puyallup$relative_importance, decreasing = TRUE),], 
          file = here("puyallup",
                      "output","puyallup_chinook_relative_importance.csv"))




# 
# data_day_night %>% 
#   filter(doy>120, doy <= 200, daytime_category == "Night") %>% 
#   ggplot() +
#   geom_line(aes(x = doy, y = resid, group = as.factor(year), col = "temperature \n residuals"),alpha = 0.5) +
#   geom_line(aes(x = doy, y = temp, group= as.factor(year), col = "temperature"),alpha = 0.5)+
#   geom_line(aes(x = doy, y = photoperiod, group = as.factor(year), col = "photoperiod"), 
#             linewidth = 2,alpha = 0.7)+
#   scale_color_manual(name = "",values = c("temperature \n residuals" = "slategray",
#                                           "temperature" = "#B17A79",
#                                           "photoperiod" = "#E8D68A")) +
#   theme_classic() +
#   labs(x = "Day of Year",
#        y = "Value",
#        title = "") +
#   theme(legend.position = "right")
# 
# ggsave(here("puyallup","output","puyallup_temperature_residuals.png"), width = 7, height = 5, units = "in", dpi = 300)  
# 

#####

# anomaly models


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


data %>% 
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
         
  
  ggplot() +
  # geom_line(aes(x = doy, y = temp_anomaly_day, col = "temperature anomaly"), alpha = 0.5) +
  # geom_line(aes(x = doy, y = temp_day, col = "temperature"), alpha = 0.5) +
  # geom_line(aes(x = doy, y = flow_anomaly_day, col = "flow anomaly"), alpha = 0.5) +
  # geom_line(aes(x = doy, y = flow_day, col = "flow"), alpha = 0.5) +
  geom_line(aes(x = doy, y = chinook0_hatchery_perhour_day_diff, col = "chinook0 hatchery"), alpha = 0.5) +
  facet_wrap(~year, scales = "free") +
  scale_color_manual(name = "", values = c("temperature anomaly" = "slategray",
                                           "flow anomaly" = "blue",
                                           "flow" = "lightblue",
                                           "temperature" = "red")) +
  theme_classic() +
  labs(x = "Day of Year",
       y = "Anomaly",
       title = "") 


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

data %>% 
  ungroup() %>%
  filter(doy > 130 & doy <= 218) %>%
  dplyr::select(flow_anomaly_day, temp_anomaly_day,
         season_day, 
         temp_diff_day, flow_diff_day) %>%
  GGally::ggpairs(aes(alpha = 0.2))

ggsave(here("puyallup","output","chinook_covariates_correlation_anomaly.png"),width = 10, height = 10)


data %>% 
  filter(doy>130, doy <= 218) %>% 
  ggplot() +
  # geom_line(aes(x = doy, y = resid, group = as.factor(year), col = "temperature \n residuals"),alpha = 0.5) +
  geom_line(aes(x = doy, y = temp_day, group= as.factor(year), col = "temperature"),alpha = 0.5)+
  geom_line(aes(x = doy, y = temp_anomaly_day, group= as.factor(year), col = "temperature \n anomaly"),alpha = 0.5)+
  # geom_line(aes(x = doy, y = photoperiod, group = as.factor(year), col = "photoperiod"), 
  # linewidth = 2,alpha = 0.7)+
  scale_color_manual(name = "",values = c("temperature \n anomaly" = "slategray",
                                          "temperature" = "#B17A79")) + 
  # "photoperiod" = "#E8D68A")) +
  theme_classic() +
  labs(x = "Day of Year",
       y = "Value",
       title = "") +
  theme(legend.position = "right")

ggsave(here("puyallup","output","puyallup_temperature_anomaly.png"), width = 7, height = 5, units = "in", dpi = 300)  

# flow anomaly

data %>% 
  filter(doy>130, doy <= 218) %>% 
  ggplot() +
  geom_line(aes(x = doy, y = flow_day, group= as.factor(year), col = "flow"),alpha = 0.5)+
  geom_line(aes(x = doy, y = flow_anomaly_day, group= as.factor(year), col = "flow \n anomaly"),alpha = 0.5)+
  scale_color_manual(name = "",values = c("flow \n anomaly" = "slategray",
                                          "flow" = "#93BCCD")) + 
  theme_classic() +
  labs(x = "Day of Year",
       y = "Value",
       title = "") +
  theme(legend.position = "right")

ggsave(here("puyallup","output","puyallup_flow_anomaly.png"), width = 7, height = 5, units = "in", dpi = 300)



covariates_chinook0_puyallup_w_temp <- arrange(data,doy) %>%
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
list_combinations <- get_covariate_combinations(1:6)
out.tab_season_puyallup <- NULL
out_riv_puyallup <- NULL
fits_season_puyallup <- list()

amount_done = length(fits_season_puyallup)
for(i in (amount_done + 1):length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  season = 0
  temperature_difference = 0
  flow_anomaly = 0
  temp_anomaly = 0
  # lunar_phase = 0
  hatchery = 0
  flow_difference = 0
  for(j in covariates){
    if(j == 1){
      k = 1
      flow_anomaly =1
    }
    else if(j==2){
      k = 2
      temp_anomaly = 1
    }
    else if(j==3){
      k = 3
      season = 1
    }
    else if(j==4){
      k = 4
      flow_difference = 1
    }
    
    else if(j==5){
      k = 5
      temperature_difference = 1
    }
    else if(j==6){
      k = 6
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_puyallup_w_temp[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_puyallup_w_temp)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  fit.model = c(list(c= c), mod_list(num_rows,c_num,
                                     ifelse(k==6,1,0), FALSE, FALSE))
  
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))


  out=data.frame(c=name, season = season,
                 temperature_difference = temperature_difference, 
                 flow_anomaly = flow_anomaly,
                 temp_anomaly = temp_anomaly,
                 flow_difference = flow_difference,
                 # lunar_phase = lunar_phase, 
                 hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv_puyallup <- rbind(out,out_riv_puyallup)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    temperature_difference = ifelse(temperature_difference==1, fit$par$U[which(covariates == 5),], NA),
                    flow_anomaly = ifelse(flow_anomaly == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_difference = ifelse(flow_difference == 1, fit$par$U[which(covariates == 4),], NA),
                    temp_anomaly = ifelse(temp_anomaly == 1, fit$par$U[which(covariates == 2),], NA),
                    # lunar_phase = ifelse(lunar_phase == 1, fit$par$U[which(covariates == 2),], NA),
                    hatchery_day = ifelse(hatchery == 1, fit$par$U["day",], NA),
                    hatchery_night = ifelse(hatchery == 1, fit$par$U["night",], NA),
                    AICc=fit$AICc)
  out.tab_season_puyallup=rbind(out.tab_season_puyallup,out2)
  fits_season_puyallup=c(fits_season_puyallup,list(fit))
  # 
  # 
}

out.tab_season_puyallup$deltaAICc <- out.tab_season_puyallup$AICc - min(out.tab_season_puyallup$AICc)
min.AICc <- order(out.tab_season_puyallup$AICc)
out.tab_season_puyallup.ordered <- out.tab_season_puyallup[min.AICc, ]
out.tab_season_puyallup.ordered

#round all the values to 2 decimal places

out.tab_season_puyallup.ordered <- round(out.tab_season_puyallup.ordered, 2)

write.csv(out.tab_season_puyallup.ordered, 
          file = here("puyallup",
                      "output","model_selection_puyallup_chinook_anomaly.csv"))


#best model

puyallup_chinook_best_model <- fits_season_puyallup[[which.min(sapply(fits_season_puyallup, function(x) x$AICc))]]

#save

save(puyallup_chinook_best_model, file = here("puyallup",
                                              "output","puyallup_chinook_best_model_anomaly.RData"))

#ci
tidy(puyallup_chinook_best_model)

autoplot(puyallup_chinook_best_model)


#relative importance

fit.model = c(mod_list(num_rows,0,0, FALSE, FALSE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temperature_difference = 0
temp_anomaly = 0
flow_anomaly = 0
# lunar_phase = 0
hatchery = 0
flow_difference = 0


out=data.frame(c=name, season = season,
               temperature_difference = temperature_difference, flow_anomaly = flow_anomaly,
               temp_anomaly = temp_anomaly, flow_difference = flow_difference,
               hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out_riv_puyallup$deltaAICc <- NULL
out_riv_puyallup <- rbind(out,out_riv_puyallup)


out_riv_puyallup$rel.LL <- NULL
out_riv_puyallup$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_season_puyallup=c(fits_season_puyallup,list(fit))

weights <- akaike.weights(out_riv_puyallup$AICc)

out_riv_puyallup$deltaAICc <- weights$deltaAIC
out_riv_puyallup$rel.LL <- weights$rel.LL
out_riv_puyallup$weights <- weights$weights


min.AICc <- order(out_riv_puyallup$AICc)
out_riv_puyallup.ordered <- out_riv_puyallup[min.AICc, ]
out_riv_puyallup.ordered

out_riv_puyallup.ordered$cumulative_weights <- cumsum(out_riv_puyallup.ordered$weights)

relative_importance_season <- sum(out_riv_puyallup$weights[out_riv_puyallup$season==1])
relative_importance_temperature_difference <- sum(out_riv_puyallup$weights[out_riv_puyallup$temperature_difference==1])
relative_importance_flow_anomaly <- sum(out_riv_puyallup$weights[out_riv_puyallup$flow_anomaly==1])
relative_importance_temp_anomaly <- sum(out_riv_puyallup$weights[out_riv_puyallup$temp_anomaly==1])
relative_importance_flow_difference <- sum(out_riv_puyallup$weights[out_riv_puyallup$flow_difference==1])
# relative_importance_lunar_phase <- sum(out_riv_puyallup$weights[out_riv_puyallup$lunar_phase==1])
relative_importance_hatchery <- sum(out_riv_puyallup$weights[out_riv_puyallup$hatchery==1])

riv_season_puyallup <- data.frame(variable = c("season",
                                               "temperature difference",
                                               "flow anomaly",
                                               "temp anomaly",
                                               "flow difference",
                                               # "lunar phase",
                                               "hatchery"),
                                  relative_importance = c(relative_importance_season,
                                                          relative_importance_temperature_difference,
                                                          relative_importance_flow_anomaly,
                                                          relative_importance_temp_anomaly,
                                                          relative_importance_flow_difference,
                                                          # relative_importance_lunar_phase,
                                                          relative_importance_hatchery))


#save
write.csv(riv_season_puyallup[order(riv_season_puyallup$relative_importance, decreasing = TRUE),], 
          file = here("puyallup",
                      "output","puyallup_chinook_relative_importance_anomaly.csv"))


# Figures

#####




predict_chinook0_puyallup <- predict(puyallup_chinook_best_model, type = "ytT", interval = "confidence")

glimpse(predict_chinook0_puyallup$pred)

head(predict_chinook0_puyallup$pred)

predict_chinook0_puyallup$pred$trap <- 'screw'

#last 4 digits of the rownames are the year
predict_chinook0_puyallup$pred$year <-  as.numeric(substr(predict_chinook0_puyallup$pred$.rownames,
                                                          nchar(predict_chinook0_puyallup$pred$.rownames)-3,
                                                          nchar(predict_chinook0_puyallup$pred$.rownames)))

predict_chinook0_puyallup$pred$daynight_category <- ifelse(substr(predict_chinook0_puyallup$pred$.rownames, 11,13) == 'day', 'day', 'night')

predict_chinook0_puyallup$pred$doy <- predict_chinook0_puyallup$pred$t+130

puyallup_covariates_chinook0 <- as.data.frame(t(covariates_chinook0_puyallup_w_temp))

puyallup_covariates_chinook0$doy <- as.numeric(rownames(puyallup_covariates_chinook0))

puyallup_covariates_chinook0_long <-  puyallup_covariates_chinook0 %>% 
  dplyr::select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","daynight_category","year"),
               names_pattern = "(.*)_(.*)_diff_(.{4})") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'day', 'day', 'night'))

predict_chinook0_puyallup$pred <- predict_chinook0_puyallup$pred %>%
  left_join(puyallup_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_chinook0_puyallup$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"), alpha = 0.8)+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2, alpha = 0.6)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_puyallup$pred, aes(x = doy, y = chinook0_hatchery_perhour, 
                                                       color = "hatchery difference"), alpha = 0.5) +
  facet_wrap(~year+daynight_category, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", 
                                           "hatchery difference" = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(130, 200), breaks = c(140,160,180))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("puyallup","output","chinook0_puyallup_prediction_w_hatchery_anomaly.png"), width = 6, height = 8, units = "in", dpi = 300)  


autoplot(puyallup_chinook_best_model, plot.type = "std.model.resids.ytT") +
  geom_point(aes(color = "wild"), size = 1)+
  labs(x = "Day of year")+
  scale_color_manual(name = "", values = c("wild" = "salmon"))+
  theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 6),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("puyallup","output","chinook0_puyallup_residuals_anomaly.png"), width = 8, height = 6, units = "in", dpi = 300)  

