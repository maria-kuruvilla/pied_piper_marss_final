# Goal - To run marss analysis on coho data with
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

data_day <- read.csv(here("data","dungeness_unaggregated_day_w_covariates.csv"))
data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_day$daytime_category <- "Day"
data_night$daytime_category <- "Night"


#interpolate the hatchery covariate values
data_day$coho1_hatchery_perhour_interpolate <- na.approx(data_day$coho1_hatchery_perhour, na.rm = FALSE)
data_night$coho1_hatchery_perhour_interpolate <- na.approx(data_night$coho1_hatchery_perhour, na.rm = FALSE)


data_day_night <- rbind(data_day,data_night)

#calculate the cumulative sum of coho1_per_hour

data_day_night <- data_day_night %>%
  filter(doy > 100, doy < 200) %>%
  group_by(year, daytime_category) %>% 
  mutate(coho1_wild_perhour_cumsum = cumsum(ifelse(is.na(coho1_wild_perhour),
                                                      0,
                                                      coho1_wild_perhour))) %>% 
  mutate(coho1_proportion = coho1_wild_perhour_cumsum/
           sum(coho1_wild_perhour,na.rm = TRUE))
#calculate the median day of migration for each year (day of year that 0.5 proportion of fish have passed)

# make column for median day of migration - day of year and call it season

coho_median <- data_day_night %>% 
  group_by(year, daytime_category) %>% 
  dplyr::select(year,coho1_proportion,doy) %>%
  summarise(coho_median_doy = doy[which.min(abs(coho1_proportion - 0.5))])

data_day_night <- left_join(data_day_night,coho_median,by = c("year","daytime_category"))

#calculate the difference between the median day of migration and the day of year

data_day_night <- data_day_night %>%
  mutate(season = coho_median_doy - doy)

data_day_night %>%
  filter(daytime_category == "Night") %>%
  ggplot() +
  geom_line(aes(x = doy, y = scale(coho1_wild_perhour, center = FALSE), group = year))+
  geom_line(aes(x = doy, y = scale(season, center = FALSE), group = year))+
  # geom_line(aes(x = doy, y = scale(lunar_phase, center = FALSE), group = year))+
  geom_line(aes(x = doy, y = scale(flow, center = TRUE), group = year))+
  scale_color_distiller(palette = "Set1")+
  facet_wrap(~year, scale = "free_y")

data_day_night %>%
  filter(daytime_category == "Day" & doy > 120 & doy <= 200) %>%
  ggplot() +
  # geom_line(aes(x = doy, y = scale(photoperiod), group = year), col = "darkred") +
  # geom_line(aes(x = doy, y = season, group = year), col = "darkblue") +
  # geom_line(aes(x = doy, y = scale(season, center = FALSE), group = year), col = "darkred") +
  # ylim(-1,1)+
  # geom_line(aes(x = doy, y = photo_diff, group = year), col = "darkblue") +
  # geom_line(aes(x = doy, y = temp_diff, group = year), col = "darkgreen") +
  # geom_line(aes(x = doy, y = flow_diff, group = year), col = "darkorange") +
  # geom_line(aes(x = doy, y = flow, group = year, col = as.factor(year)), linewidth  = 1, alpha = 0.4) +
  geom_line(aes(x = doy, y = temp, group = year, col = as.factor(year)), linewidth  = 1, alpha = 0.4) +
  # geom_line(aes(x = doy, y = lunar_phase, group = year), col = "darkgrey") +
  # geom_line(aes(x = doy, y = resid, group = year), col = "black") +
  theme_minimal()

data_day_night %>% 
  ungroup() %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  select(flow,lunar_phase, season, resid, temp_diff, flow_diff, photo_diff) %>%
  GGally::ggpairs(aes(alpha = 0.2))

ggsave(here("dungeness","output","coho_covariates_correlation.png"),width = 10, height = 10)

covariates_coho1 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, daytime_category, flow,
                lunar_phase, season, 
                resid, temp_diff, flow_diff, photo_diff, 
                coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), 
              values_from = c(flow, lunar_phase, season, 
                              resid, temp_diff, flow_diff, photo_diff, 
                              coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2
total_covariates = dim(covariates_coho1)[1]
covariate_num = total_covariates/num_rows

#scale and center flow and lunar phase
for(i in (1):(num_rows*2)){
  if(sum(covariates_coho1[i,]) != 0){
    print(rownames(covariates_coho1)[i])
    covariates_coho1[i,] = scale(covariates_coho1[i,],
                                    center = TRUE, scale= TRUE)[,1]
  }
}

for(i in (num_rows*2 + 1):(total_covariates)){
  if(sum(covariates_coho1[i,]) != 0){
    print(rownames(covariates_coho1)[i])
    covariates_coho1[i,] = scale(covariates_coho1[i,],
                                    center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_coho_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}


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


#######
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
fit_equal_q_coho <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                     control=list(maxit=2000))
fit_equal_q_coho$AICc

fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,TRUE))
fit_unequal_q_coho <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                       control=list(maxit=2000))
fit_unequal_q_coho$AICc

df_errors_coho <- data.frame(Error_structure = c("Equal","Unequal"), 
                        AICc = c(fit_equal_q_coho$AICc,fit_unequal_q_coho$AICc))

#save table
write.csv(df_errors_coho, file = here("dungeness",
                                 "output","error_structure_dungeness_coho.csv"))

#correalted covariates - photo_diff and season
out.tab.all.years.equalq <- NULL
fits.all.years.equalq <- NULL
nyears = num_years*2
for(kk in c(3,7)){
  c = covariates_coho1[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0,FALSE,FALSE))
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  ci = tidy(fit)
  out=data.frame(c=name_individual, Estimate = ci[33,2],
                 conf_low = ci[33,4], conf_high = ci[33,5],
                 AICc=fit$AICc)
  out.tab.all.years.equalq=rbind(out.tab.all.years.equalq,out)
  fits.all.years.equalq=c(fits.all.years.equalq,list(fit))
}

out.tab.all.years.equalq$delta_AICc <- out.tab.all.years.equalq$AICc  - min(out.tab.all.years.equalq$AICc)

out.tab.all.years.equalq <- out.tab.all.years.equalq[order(out.tab.all.years.equalq$AICc),]

write.csv(out.tab.all.years.equalq, file = here("dungeness",
                                                  "output","correlated_covariates_dungenss_coho.csv"))






num_years = 15 
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:7)
out.tab_season_coho <- NULL
out_riv_coho <- NULL
fits_season_coho <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  season = 0
  resid = 0
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
      k = 2
      lunar_phase = 1
    }
    else if(j==3){
      k = 3
      season = 1
    }
    else if(j==4){
      k = 4
      resid = 1
    }
    else if(j==5){
      k = 5
      temperature_difference = 1
    }
    
    else if(j==6){
      k = 6
      flow_difference = 1
    }
    else if(j==7){
      k = 8
      hatchery = 1
    }
    
    c = rbind(c,covariates_coho1[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_coho1)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==8){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, FALSE))
  } else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, FALSE))
  }
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season, resid = resid,
                 temperature_difference = temperature_difference, flow = flow,
                 flow_difference = flow_difference,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv_coho <- rbind(out,out_riv_coho)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    resid = ifelse(resid==1, fit$par$U[which(covariates == 4),], NA),
                    temperature_difference = ifelse(temperature_difference==1, fit$par$U[which(covariates == 5),], NA),
                    flow = ifelse(flow == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_difference = ifelse(flow_difference == 1, fit$par$U[which(covariates == 6),], NA),
                    lunar_phase = ifelse(lunar_phase == 1, fit$par$U[which(covariates == 2),], NA),
                    hatchery_day = ifelse(hatchery == 1, fit$par$U["day",], NA),
                    hatchery_night = ifelse(hatchery == 1, fit$par$U["night",], NA),
                    AICc=fit$AICc)
  out.tab_season_coho=rbind(out.tab_season_coho,out2)
  fits_season_coho=c(fits_season_coho,list(fit))
  
  
}

out.tab_season_coho$deltaAICc <- out.tab_season_coho$AICc - min(out.tab_season_coho$AICc)
min.AICc <- order(out.tab_season_coho$AICc)
out.tab_season_coho.ordered <- out.tab_season_coho[min.AICc, ]
out.tab_season_coho.ordered

write.csv(out.tab_season_coho.ordered, 
          file = here("dungeness",
                      "output","model_selection_dungenss_coho.csv"))



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
resid = 0


out=data.frame(c=name, season = season, resid = resid,
               temperature_difference = temperature_difference, flow = flow,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out_riv_coho <- rbind(out,out_riv_coho)

# out_riv$deltaAICc <- NULL
out_riv_coho$rel.LL <- NULL
out_riv_coho$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_season_coho=c(fits_season_coho,list(fit))

weights <- akaike.weights(out_riv_coho$AICc)

out_riv_coho$deltaAICc <- weights$deltaAIC
out_riv_coho$rel.LL <- weights$rel.LL
out_riv_coho$weights <- weights$weights


min.AICc <- order(out_riv_coho$AICc)
out_riv_coho.ordered <- out_riv_coho[min.AICc, ]
out_riv_coho.ordered

out_riv_coho.ordered$cumulative_weights <- cumsum(out_riv_coho.ordered$weights)

relative_importance_season <- sum(out_riv_coho$weights[out_riv_coho$season==1])
realative_importance_resid <- sum(out_riv_coho$weights[out_riv_coho$resid==1])
relative_importance_temperature_difference <- sum(out_riv_coho$weights[out_riv_coho$temperature_difference==1])
relative_importance_flow <- sum(out_riv_coho$weights[out_riv_coho$flow==1])
relative_importance_flow_difference <- sum(out_riv_coho$weights[out_riv_coho$flow_difference==1])
relative_importance_lunar_phase <- sum(out_riv_coho$weights[out_riv_coho$lunar_phase==1])
relative_importance_hatchery <- sum(out_rivv$weights[out_riv_coho$hatchery==1])

riv_season_coho <- data.frame(variable = c("season",
                                           "resid",
                                      "temperature difference",
                                      "flow",
                                      "flow difference",
                                      "lunar phase",
                                      "hatchery"),
                         relative_importance = c(relative_importance_season,
                                                 realative_importance_resid,
                                                 relative_importance_temperature_difference,
                                                 relative_importance_flow,
                                                 relative_importance_flow_difference,
                                                 relative_importance_lunar_phase,
                                                 relative_importance_hatchery))


#save
write.csv(riv_season_coho, 
          file = here("dungeness",
                      "output","dungenss_coho_relative_importance.csv"))

#best model

dungeness_coho_best_model <- fits_season_coho[[which.min(sapply(fits_season_coho, function(x) x$AICc))]]

#save

save(dungeness_coho_best_model, file = here("dungeness",
                                               "output","dungeness_coho_best_model.RData"))

#ci
tidy(dungeness_coho_best_model)

autoplot(dungeness_coho_best_model)

# Figures

#####



predict_coho1_dungeness <- predict(dungeness_coho_best_model, type = "ytT", interval = "confidence")

glimpse(predict_coho1_dungeness$pred)

head(predict_coho1_dungeness$pred)

predict_coho1_dungeness$pred$trap <- 'screw'

predict_coho1_dungeness$pred$year <-  as.numeric(substr(predict_coho1_dungeness$pred$.rownames, 
                                                           1, 4))

predict_coho1_dungeness$pred$daynight_category <- ifelse(substr(predict_coho1_dungeness$pred$.rownames, 
                                                                   6, nchar(predict_coho1_dungeness$pred$.rownames)) == 'Day', 'day', 'night')

predict_coho1_dungeness$pred$doy <- predict_coho1_dungeness$pred$t+120

dungeness_covariates_coho1 <- as.data.frame(t(covariates_coho1))

dungeness_covariates_coho1$doy <- as.numeric(rownames(dungeness_covariates_coho1))

dungeness_covariates_coho1_long <-  dungeness_covariates_coho1 %>% 
  select(doy, starts_with("coho")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_coho1_dungeness$pred <- predict_coho1_dungeness$pred %>%
  left_join(dungeness_covariates_coho1_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_coho1_dungeness$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_dungeness$pred, aes(x = doy, y = log(coho1_hatchery_perhour_interpolate+1), 
                                                        color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 5, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
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

ggsave(here("dungeness","output","coho1_dungeness_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)  


autoplot(dungeness_coho_best_model, plot.type = "std.model.resids.ytT") +
  geom_point(aes(color = "wild"), size = 1.5)+
  labs(x = "Day of year")+
  scale_color_manual(name = "", values = c("wild" = "salmon"))+
  theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("dungeness","output","coho1_dungeness_residuals.png"), width = 8, height = 6, units = "in", dpi = 300)  


#####


#just do analysis for night values

data_night <- read.csv(here("data","dungeness_unaggregated_night_w_covariates.csv"))

data_night$daytime_category = "Night"
data_night$coho1_hatchery_perhour_interpolate <- na.approx(data_night$coho1_hatchery_perhour, na.rm = FALSE)

#calculate the cumulative sum of coho1_per_hour

data_night <- data_night %>%
  filter(doy > 100, doy < 200) %>%
  group_by(year) %>% 
  mutate(coho1_wild_perhour_cumsum = cumsum(ifelse(is.na(coho1_wild_perhour),
                                                   0,
                                                   coho1_wild_perhour))) %>% 
  mutate(coho1_proportion = coho1_wild_perhour_cumsum/
           sum(coho1_wild_perhour,na.rm = TRUE))
#calculate the median day of migration for each year (day of year that 0.5 proportion of fish have passed)

# make column for median day of migration - day of year and call it season

coho_median_night <- data_night %>% 
  group_by(year) %>% 
  select(year,coho1_proportion,doy) %>%
  summarise(coho_median_doy = doy[which.min(abs(coho1_proportion - 0.5))])

data_night <- left_join(data_night,coho_median_night,by = c("year"))

#calculate the difference between the median day of migration and the day of year

data_night <- data_night %>%
  mutate(season = coho_median_doy - doy)

data_night %>%
  # filter(daytime_category == "Night") %>%
  ggplot() +
  geom_line(aes(x = doy, y = scale(coho1_wild_perhour, center = FALSE), group = year))+
  geom_line(aes(x = doy, y = scale(season, center = FALSE), group = year))+
  # geom_line(aes(x = doy, y = scale(lunar_phase, center = FALSE), group = year))+
  geom_line(aes(x = doy, y = scale(flow, center = TRUE), group = year))+
  scale_color_distiller(palette = "Set1")+
  facet_wrap(~year, scale = "free_y")


covariates_coho1_night <- arrange(data_night,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, daytime_category, flow,
                lunar_phase, season, 
                resid, temp_diff, flow_diff, photo_diff, 
                coho1_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), 
              values_from = c(flow, lunar_phase, season, 
                              resid, temp_diff, flow_diff, photo_diff, 
                              coho1_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years
total_covariates = dim(covariates_coho1_night)[1]
covariate_num = total_covariates/num_rows

#scale and center flow and lunar phase
for(i in (1):(num_rows*2)){
  if(sum(covariates_coho1_night[i,]) != 0){
    print(rownames(covariates_coho1_night)[i])
    covariates_coho1_night[i,] = scale(covariates_coho1_night[i,],
                                 center = TRUE, scale= TRUE)[,1]
  }
}

for(i in (num_rows*2 + 1):(total_covariates)){
  if(sum(covariates_coho1_night[i,]) != 0){
    print(rownames(covariates_coho1_night)[i])
    covariates_coho1_night[i,] = scale(covariates_coho1_night[i,],
                                 center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_coho_summer_perhour <- arrange(data_night,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}



Cmat_skagit_night <- function(nyears,ncov,hatchery=0){
  vars = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")
  
  
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
  
  return(C)
}


#######





#function for mod list

######

mod_list <- function(nyears,ncov,hatchery=0){
  
  Q = "diagonal and equal"
  
  
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
      C = Cmat_skagit_night(nyears,ncov,hatchery)
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
#####



#correalted covariates - photo_diff and season
out.tab.all.years.equalq <- NULL
fits.all.years.equalq <- NULL
nyears = num_years
for(kk in c(3,7)){
  c = covariates_coho1_night[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_coho1_night)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-11)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0))
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  ci = tidy(fit)
  out=data.frame(c=name_individual, Estimate = ci[18,2],
                 conf_low = ci[18,4], conf_high = ci[18,5],
                 AICc=fit$AICc)
  out.tab.all.years.equalq=rbind(out.tab.all.years.equalq,out)
  fits.all.years.equalq=c(fits.all.years.equalq,list(fit))
}

out.tab.all.years.equalq$delta_AICc <- out.tab.all.years.equalq$AICc  - min(out.tab.all.years.equalq$AICc)

out.tab.all.years.equalq <- out.tab.all.years.equalq[order(out.tab.all.years.equalq$AICc),]

write.csv(out.tab.all.years.equalq, file = here("dungeness",
                                                "output","correlated_covariates_dungenss_coho_night.csv"))



num_years = 15 
num_rows = num_years
list_combinations <- get_covariate_combinations(1:7)
out.tab_season_coho <- NULL
out_riv_coho <- NULL
fits_season_coho <- list()

amount_done <- length(fits_season_coho)
for(i in (amount_done + 1):length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  season = 0
  resid = 0
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
      k = 2
      lunar_phase = 1
    }
    else if(j==3){
      k = 3
      season = 1
    }
    else if(j==4){
      k = 4
      resid = 1
    }
    else if(j==5){
      k = 5
      temperature_difference = 1
    }
    
    else if(j==6){
      k = 6
      flow_difference = 1
    }
    else if(j==7){
      k = 8
      hatchery = 1
    }
    
    c = rbind(c,covariates_coho1_night[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_coho1_night)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-11))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==8){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  } else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  }
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season, resid = resid,
                 temperature_difference = temperature_difference, flow = flow,
                 flow_difference = flow_difference,
                 lunar_phase = lunar_phase, hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv_coho <- rbind(out,out_riv_coho)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    resid = ifelse(resid==1, fit$par$U[which(covariates == 4),], NA),
                    temperature_difference = ifelse(temperature_difference==1, fit$par$U[which(covariates == 5),], NA),
                    flow = ifelse(flow == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_difference = ifelse(flow_difference == 1, fit$par$U[which(covariates == 6),], NA),
                    lunar_phase = ifelse(lunar_phase == 1, fit$par$U[which(covariates == 2),], NA),
                    hatchery = ifelse(hatchery == 1, fit$par$U[which(covariates == 7),], NA),
                    AICc=fit$AICc)
  out.tab_season_coho=rbind(out.tab_season_coho,out2)
  fits_season_coho=c(fits_season_coho,list(fit))
  
  
}

out.tab_season_coho$deltaAICc <- out.tab_season_coho$AICc - min(out.tab_season_coho$AICc)
min.AICc <- order(out.tab_season_coho$AICc)
out.tab_season_coho.ordered <- out.tab_season_coho[min.AICc, ]
out.tab_season_coho.ordered <- round(out.tab_season_coho.ordered, 2)

write.csv(out.tab_season_coho.ordered, 
          file = here("dungeness",
                      "output","model_selection_dungenss_coho_night.csv"))


#best model

dungeness_coho_best_model_night <- fits_season_coho[[which.min(sapply(fits_season_coho, function(x) x$AICc))]]

#save

save(dungeness_coho_best_model_night, file = here("dungeness",
                                            "output","dungeness_coho_best_model_night.RData"))

#ci
tidy(dungeness_coho_best_model_night)

autoplot(dungeness_coho_best_model_night)

#relative importance

fit.model = c(mod_list(num_rows,0,0))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
resid = 0


out=data.frame(c=name, season = season, resid = resid,
               temperature_difference = temperature_difference, flow = flow,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out_riv_coho <- rbind(out,out_riv_coho)

# out_riv$deltaAICc <- NULL
out_riv_coho$rel.LL <- NULL
out_riv_coho$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_season_coho=c(fits_season_coho,list(fit))

weights <- akaike.weights(out_riv_coho$AICc)

out_riv_coho$deltaAICc <- weights$deltaAIC
out_riv_coho$rel.LL <- weights$rel.LL
out_riv_coho$weights <- weights$weights


min.AICc <- order(out_riv_coho$AICc)
out_riv_coho.ordered <- out_riv_coho[min.AICc, ]
out_riv_coho.ordered

out_riv_coho.ordered$cumulative_weights <- cumsum(out_riv_coho.ordered$weights)

relative_importance_season <- sum(out_riv_coho$weights[out_riv_coho$season==1])
realative_importance_resid <- sum(out_riv_coho$weights[out_riv_coho$resid==1])
relative_importance_temperature_difference <- sum(out_riv_coho$weights[out_riv_coho$temperature_difference==1])
relative_importance_flow <- sum(out_riv_coho$weights[out_riv_coho$flow==1])
relative_importance_flow_difference <- sum(out_riv_coho$weights[out_riv_coho$flow_difference==1])
relative_importance_lunar_phase <- sum(out_riv_coho$weights[out_riv_coho$lunar_phase==1])
relative_importance_hatchery <- sum(out_riv_coho$weights[out_riv_coho$hatchery==1])

riv_season_coho <- data.frame(variable = c("season",
                                           "resid",
                                           "temperature difference",
                                           "flow",
                                           "flow difference",
                                           "lunar phase",
                                           "hatchery"),
                              relative_importance = c(relative_importance_season,
                                                      realative_importance_resid,
                                                      relative_importance_temperature_difference,
                                                      relative_importance_flow,
                                                      relative_importance_flow_difference,
                                                      relative_importance_lunar_phase,
                                                      relative_importance_hatchery))


#save
write.csv(riv_season_coho, 
          file = here("dungeness",
                      "output","dungenss_coho_relative_importance_night.csv"))


# FIGURES





predict_coho1_dungeness_night <- predict(dungeness_coho_best_model_night, type = "ytT", interval = "confidence")

glimpse(predict_coho1_dungeness_night$pred)


head(predict_coho1_dungeness_night$pred)

predict_coho1_dungeness_night$pred$trap <- 'screw'

predict_coho1_dungeness_night$pred$year <-  as.numeric(substr(predict_coho1_dungeness_night$pred$.rownames, 
                                                        1, 4))

predict_coho1_dungeness_night$pred$daynight_category <- ifelse(substr(predict_coho1_dungeness_night$pred$.rownames, 
                                                                6, nchar(predict_coho1_dungeness_night$pred$.rownames)) == 'Day', 'day', 'night')

predict_coho1_dungeness_night$pred$doy <- predict_coho1_dungeness_night$pred$t+120

dungeness_covariates_coho1 <- as.data.frame(t(covariates_coho1_night))

dungeness_covariates_coho1$doy <- as.numeric(rownames(dungeness_covariates_coho1))

dungeness_covariates_coho1_long <-  dungeness_covariates_coho1 %>% 
  select(doy, starts_with("coho")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_coho1_dungeness_night$pred <- predict_coho1_dungeness_night$pred %>%
  left_join(dungeness_covariates_coho1_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_coho1_dungeness_night$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_dungeness_night$pred, aes(x = doy, y = log(coho1_hatchery_perhour_interpolate+1), 
                                                     color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 3, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(120, 160), breaks = c(130,150))+
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

ggsave(here("dungeness","output","coho1_dungeness_prediction_w_hatchery_night.png"), width = 6, height = 8, units = "in", dpi = 300)  


autoplot(dungeness_coho_best_model_night, plot.type = "std.model.resids.ytT") +
  geom_point(aes(color = "wild"), size = 1)+
  labs(x = "Day of year")+
  scale_color_manual(name = "", values = c("wild" = "salmon"))+
  theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("dungeness","output","coho1_dungeness_residuals_night.png"), width = 8, height = 6, units = "in", dpi = 300)  



# anomaly models

data_day_night_anomaly <- read.csv(here("data","dungeness_anomaly.csv"))

glimpse(data_day_night_anomaly)

data_day_night_anomaly <- data_day_night_anomaly %>%
  filter(doy > 100, doy < 200) %>%
  group_by(year, daytime_category) %>% 
  mutate(coho1_wild_perhour_cumsum = cumsum(ifelse(is.na(coho1_wild_perhour),
                                                   0,
                                                   coho1_wild_perhour))) %>% 
  mutate(coho1_proportion = coho1_wild_perhour_cumsum/
           sum(coho1_wild_perhour,na.rm = TRUE))
#calculate the median day of migration for each year (day of year that 0.5 proportion of fish have passed)

# make column for median day of migration - day of year and call it season

coho_median <- data_day_night_anomaly %>% 
  group_by(year, daytime_category) %>% 
  dplyr::select(year,coho1_proportion,doy) %>%
  summarise(coho_median_doy = doy[which.min(abs(coho1_proportion - 0.5))])

data_day_night_anomaly <- left_join(data_day_night_anomaly,coho_median,by = c("year","daytime_category"))

#calculate the difference between the median day of migration and the day of year

data_day_night_anomaly <- data_day_night_anomaly %>%
  mutate(season = coho_median_doy - doy)

covariates_coho1_anomaly <- arrange(data_day_night_anomaly,doy) %>%
  filter(year != 2015 & year != 2016 & doy >120 & doy <= 160) %>%
  dplyr::select(year,doy, daytime_category, flow_anomaly,
                temp_anomaly, season, 
                coho1_hatchery_perhour_diff) %>%
  pivot_wider(names_from = c(year, daytime_category), 
              values_from = c(flow_anomaly, temp_anomaly, season, 
                              coho1_hatchery_perhour_diff)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005-1
num_rows = num_years*2
total_covariates = dim(covariates_coho1_anomaly)[1]
covariate_num = total_covariates/num_rows

#scale and center flow and lunar phase
for(i in (1):(num_rows*2)){
  if(sum(covariates_coho1_anomaly[i,]) != 0){
    print(rownames(covariates_coho1_anomaly)[i])
    covariates_coho1_anomaly[i,] = scale(covariates_coho1_anomaly[i,],
                                 center = TRUE, scale= TRUE)[,1]
  }
}

for(i in (num_rows*2 + 1):(total_covariates)){
  if(sum(covariates_coho1_anomaly[i,]) != 0){
    print(rownames(covariates_coho1_anomaly)[i])
    covariates_coho1_anomaly[i,] = scale(covariates_coho1_anomaly[i,],
                                 center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_coho_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & year != 2016 & doy > 120 & doy <= 160) %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour)[1]){
  subset_coho_summer_perhour[i,] = scale(subset_coho_summer_perhour[i,])[,1]
}


nyears = num_years*2
c = NULL
fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,FALSE))
fit_equal_q_coho <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                          control=list(maxit=2000))
fit_equal_q_coho$AICc

fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,TRUE))
fit_unequal_q_coho <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
                            control=list(maxit=2000))
fit_unequal_q_coho$AICc

df_errors_coho <- data.frame(Error_structure = c("Equal","Unequal"), 
                             AICc = c(fit_equal_q_coho$AICc,fit_unequal_q_coho$AICc))


data_day_night_anomaly %>% 
  ungroup() %>%
  filter(year != 2015 & year != 2016 & doy > 120 & doy <= 160) %>%
  dplyr::select(season,  temp_anomaly, flow_anomaly) %>%
  GGally::ggpairs(aes(alpha = 0.2))

ggsave(here("dungeness","output","coho_covariates_anomaly_correlation.png"),width = 10, height = 10)




num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:4)
out.tab_anomaly_coho <- NULL
out_riv_coho <- NULL
fits_anomaly_coho <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  season = 0
  temp_anomaly = 0
  flow_anomaly = 0
  hatchery = 0
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
      hatchery = 1
    }
    
    
    c = rbind(c,covariates_coho1_anomaly[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_coho1_anomaly)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==4){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, FALSE))
  } else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, FALSE))
  }
  fit <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season, 
                 temp_anomaly = temp_anomaly, flow_anomaly = flow_anomaly,
                 hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv_coho <- rbind(out,out_riv_coho)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    temp_anomaly = ifelse(temp_anomaly==1, fit$par$U[which(covariates == 2),], NA),
                    flow_anomaly = ifelse(flow_anomaly == 1, fit$par$U[which(covariates == 1),], NA),
                    hatchery_day = ifelse(hatchery == 1, fit$par$U["day",], NA),
                    hatchery_night = ifelse(hatchery == 1, fit$par$U["night",], NA),
                    AICc=fit$AICc)
  out.tab_anomaly_coho=rbind(out.tab_anomaly_coho,out2)
  fits_anomaly_coho=c(fits_anomaly_coho,list(fit))
  
  
}

out.tab_anomaly_coho$deltaAICc <- out.tab_anomaly_coho$AICc - min(out.tab_anomaly_coho$AICc)
min.AICc <- order(out.tab_anomaly_coho$AICc)
out.tab_anomaly_coho.ordered <- out.tab_anomaly_coho[min.AICc, ]
out.tab_anomaly_coho.ordered

write.csv(out.tab_anomaly_coho.ordered, 
          file = here("dungeness",
                      "output","model_selection_dungenss_coho_anomaly.csv"))

#best model

dungeness_coho_best_model_anomaly <- fits_anomaly_coho[[which.min(sapply(fits_anomaly_coho, function(x) x$AICc))]]

#save

save(dungeness_coho_best_model_anomaly, file = here("dungeness",
                                                  "output","dungeness_coho_best_model_anomaly.RData"))

#ci
tidy(dungeness_coho_best_model_anomaly)

autoplot(dungeness_coho_best_model_anomaly)

# residuals for day have pattern. need to do only with night values

#relative importance

fit.model = c(mod_list(num_rows,0,0, FALSE, FALSE))
fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temp_anomaly = 0
flow_anomaly = 0
hatchery = 0



out=data.frame(c=name, season = season, 
               temp_anomaly = temp_anomaly, flow_anomaly = flow_anomaly,
               hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out_riv_coho <- rbind(out,out_riv_coho)

# out_riv$deltaAICc <- NULL
out_riv_coho$rel.LL <- NULL
out_riv_coho$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_anomaly_coho=c(fits_anomaly_coho,list(fit))

weights <- akaike.weights(out_riv_coho$AICc)

out_riv_coho$deltaAICc <- weights$deltaAIC
out_riv_coho$rel.LL <- weights$rel.LL
out_riv_coho$weights <- weights$weights


min.AICc <- order(out_riv_coho$AICc)
out_riv_coho.ordered <- out_riv_coho[min.AICc, ]
out_riv_coho.ordered

out_riv_coho.ordered$cumulative_weights <- cumsum(out_riv_coho.ordered$weights)

relative_importance_season <- sum(out_riv_coho$weights[out_riv_coho$season==1])
relative_importance_temp_anomaly <- sum(out_riv_coho$weights[out_riv_coho$temp_anomaly==1])
relative_importance_flow_anomaly <- sum(out_riv_coho$weights[out_riv_coho$flow_anomaly==1])
relative_importance_hatchery <- sum(out_riv_coho$weights[out_riv_coho$hatchery==1])

riv_anomaly_coho <- data.frame(variable = c("season",
                                           "temperature anomaly",
                                           "flow anomaly",
                                           "hatchery"),
                              relative_importance = c(relative_importance_season,
                                                      relative_importance_temp_anomaly,
                                                      relative_importance_flow_anomaly,
                                                      relative_importance_hatchery))
#round to 2 decimal places

riv_anomaly_coho$relative_importance <- round(riv_anomaly_coho$relative_importance, 2)

# order in decreasing order

riv_anomaly_coho <- riv_anomaly_coho[order(riv_anomaly_coho$relative_importance, decreasing = TRUE),]


#save
write.csv(riv_anomaly_coho, 
          file = here("dungeness",
                      "output","dungenss_coho_relative_importance_anomaly.csv"))


# FIGURES





predict_coho1_dungeness_anomaly <- predict(dungeness_coho_best_model_anomaly, type = "ytT", interval = "confidence")

glimpse(predict_coho1_dungeness_anomaly$pred)


head(predict_coho1_dungeness_anomaly$pred)

predict_coho1_dungeness_anomaly$pred$trap <- 'screw'

predict_coho1_dungeness_anomaly$pred$year <-  as.numeric(substr(predict_coho1_dungeness_anomaly$pred$.rownames, 
                                                              1, 4))

predict_coho1_dungeness_anomaly$pred$daynight_category <- ifelse(substr(predict_coho1_dungeness_anomaly$pred$.rownames, 
                                                                      6, nchar(predict_coho1_dungeness_anomaly$pred$.rownames)) == 'Day', 'day', 'night')

predict_coho1_dungeness_anomaly$pred$doy <- predict_coho1_dungeness_anomaly$pred$t+120

dungeness_covariates_coho1 <- as.data.frame(t(covariates_coho1_anomaly))

dungeness_covariates_coho1$doy <- as.numeric(rownames(dungeness_covariates_coho1))

dungeness_covariates_coho1_long <-  dungeness_covariates_coho1 %>% 
  dplyr::select(doy, starts_with("coho")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_coho1_dungeness_anomaly$pred <- predict_coho1_dungeness_anomaly$pred %>%
  left_join(dungeness_covariates_coho1_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_coho1_dungeness_anomaly$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"), alpha =0.8)+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2, alpha = 0.6)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_dungeness_anomaly$pred, aes(x = doy, y = coho1_hatchery_perhour_diff, 
                                                           color = "hatchery"), alpha =0.5) +
  facet_wrap(~year+daynight_category, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(120, 160), breaks = c(130,150))+
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

ggsave(here("dungeness","output","coho1_dungeness_prediction_w_hatchery_anomaly.png"), width = 6, height = 8, units = "in", dpi = 300)  


autoplot(dungeness_coho_best_model_anomaly, plot.type = "std.model.resids.ytT") +
  geom_point(aes(color = "wild"), size = 1)+
  labs(x = "Day of year")+
  scale_color_manual(name = "", values = c("wild" = "salmon"))+
  theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("dungeness","output","coho1_dungeness_residuals_anomaly.png"), width = 8, height = 6, units = "in", dpi = 300)  



# only night values


covariates_coho1_anomaly_night <- arrange(data_day_night_anomaly,doy) %>%
  filter(year != 2015 & doy >120 & doy <= 160, daytime_category == "Night") %>%
  dplyr::select(year,doy, daytime_category, flow_anomaly,
                temp_anomaly, season, flow_diff, temp_diff,
                coho1_hatchery_perhour_diff) %>%
  pivot_wider(names_from = c(year, daytime_category), 
              values_from = c(flow_anomaly, temp_anomaly, season, flow_diff, temp_diff,
                              coho1_hatchery_perhour_diff)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years
total_covariates = dim(covariates_coho1_anomaly_night)[1]
covariate_num = total_covariates/num_rows

#scale and center flow and lunar phase
for(i in (1):(num_rows*2)){
  if(sum(covariates_coho1_anomaly_night[i,]) != 0){
    print(rownames(covariates_coho1_anomaly_night)[i])
    covariates_coho1_anomaly_night[i,] = scale(covariates_coho1_anomaly_night[i,],
                                         center = TRUE, scale= TRUE)[,1]
  }
}

for(i in (num_rows*2 + 1):(total_covariates)){
  if(sum(covariates_coho1_anomaly_night[i,]) != 0){
    print(rownames(covariates_coho1_anomaly_night)[i])
    covariates_coho1_anomaly_night[i,] = scale(covariates_coho1_anomaly_night[i,],
                                         center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_coho_summer_perhour_night <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 120 & doy <= 160, daytime_category == "Night") %>%
  mutate(log.value = log(coho1_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
  subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
}

# no need for this because it will be equal errors

# nyears = num_years
# c = NULL
# fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,FALSE))
# fit_equal_q_coho <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
#                           control=list(maxit=2000))
# fit_equal_q_coho$AICc
# 
# fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,TRUE))
# fit_unequal_q_coho <- MARSS(subset_coho_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
#                             control=list(maxit=2000))
# fit_unequal_q_coho$AICc
# 
# df_errors_coho <- data.frame(Error_structure = c("Equal","Unequal"), 
#                              AICc = c(fit_equal_q_coho$AICc,fit_unequal_q_coho$AICc))

data_day_night_anomaly %>% 
  ungroup() %>%
  filter(year != 2015  & doy > 120 & doy <= 160, daytime_category == "Night") %>%
  dplyr::select(season,  temp_anomaly, flow_anomaly, flow_diff, temp_diff) %>%
  GGally::ggpairs(aes(alpha = 0.2))+
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)

ggsave(here("dungeness","output","coho_covariates_anomaly_night_correlation.png"),width = 10, height = 10)

# 
# 
# out.tab <- NULL
# fits <- NULL
# nyears = num_years*2
# for(kk in c(1,2,3,5)){
#   c = covariates_coho1_anomaly_night[((1+(kk-1)*nyears):(kk*nyears)),]
#   name_long = rownames(covariates_coho1_anomaly_night)[1+(kk-1)*nyears]
#   name_individual = substr(name_long,1,nchar(name_long)-9)
#   print(name_individual)
#   fit.model = c(list(c= c), mod_list(nyears,1,0,FALSE,TRUE))
#   fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
#                control=list(maxit=2000))
#   ci = tidy(fit)
#   out=data.frame(c=name_individual, Estimate = ci[34,2],
#                  conf_low = ci[34,4], conf_high = ci[34,5],
#                  AICc=fit$AICc)
#   out.tab=rbind(out.tab,out)
#   fits=c(fits,list(fit))
# }
# 
# out.tab$delta_AICc <- out.tab$AICc  - min(out.tab$AICc)
# 
# out.tab <- out.tab[order(out.tab$AICc),]
# 
# out.tab
# write.csv(out.tab, file = here("dungeness","output",
#                                "correlated_covariates_dungenss_coho_anomaly_night.csv"))


num_rows = num_years
list_combinations <- get_covariate_combinations(1:6)
out.tab_anomaly_night_coho <- NULL
out_riv_coho <- NULL
fits_anomaly_night_coho <- list()
amount_done <- length(fits_anomaly_night_coho)
for(i in (amount_done + 1):length(list_combinations)){
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  season = 0
  temp_anomaly = 0
  flow_anomaly = 0
  flow_diff = 0
  temp_diff = 0
  hatchery = 0
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
      flow_diff = 1
    }
    else if(j==5){
      k = 5
      temp_diff = 1
    }
    else if(j==6){
      k = 6
      hatchery = 1
    }
    
    
    c = rbind(c,covariates_coho1_anomaly_night[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_coho1_anomaly_night)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==4){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  } else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  }
  fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season, 
                 temp_anomaly = temp_anomaly, flow_anomaly = flow_anomaly,
                 flow_diff = flow_diff, temp_diff = temp_diff,
                 hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv_coho <- rbind(out,out_riv_coho)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    temp_anomaly = ifelse(temp_anomaly==1, fit$par$U[which(covariates == 2),], NA),
                    flow_anomaly = ifelse(flow_anomaly == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_diff = ifelse(flow_diff == 1, fit$par$U[which(covariates == 4),], NA),
                    temp_diff = ifelse(temp_diff == 1, fit$par$U[which(covariates == 5),], NA),
                    hatchery = ifelse(hatchery == 1, fit$par$U[which(covariates == 6),], NA),
                    AICc=fit$AICc)
  out.tab_anomaly_night_coho=rbind(out.tab_anomaly_night_coho,out2)
  fits_anomaly_night_coho=c(fits_anomaly_night_coho,list(fit))
  
  
}

out.tab_anomaly_night_coho$deltaAICc <- out.tab_anomaly_night_coho$AICc - min(out.tab_anomaly_night_coho$AICc)
min.AICc <- order(out.tab_anomaly_night_coho$AICc)
out.tab_anomaly_night_coho.ordered <- out.tab_anomaly_night_coho[min.AICc, ]
out.tab_anomaly_night_coho.ordered

#round the estimates to 2 decimal places

out.tab_anomaly_night_coho.ordered<- round(out.tab_anomaly_night_coho.ordered, 2)

write.csv(out.tab_anomaly_night_coho.ordered, 
          file = here("dungeness",
                      "output","model_selection_dungenss_coho_anomaly_night.csv"))

#best model

dungeness_coho_best_model_anomaly_night <- fits_anomaly_night_coho[[which.min(sapply(fits_anomaly_night_coho, function(x) x$AICc))]]

#save

save(dungeness_coho_best_model_anomaly_night, file = here("dungeness",
                                                    "output","dungeness_coho_best_model_anomaly_night.RData"))
######### start here
#ci
tidy(dungeness_coho_best_model_anomaly_night)

autoplot(dungeness_coho_best_model_anomaly_night)



#relative importance

fit.model = c(mod_list(num_rows,0,0))
fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temp_anomaly = 0
flow_anomaly = 0
hatchery = 0
temp_diff = 0
flow_diff = 0



out=data.frame(c=name, season = season, 
               temp_anomaly = temp_anomaly, flow_anomaly = flow_anomaly,
               flow_diff = flow_diff, temp_diff = temp_diff,
               hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)
out_riv_coho <- rbind(out,out_riv_coho)

# out_riv$deltaAICc <- NULL
out_riv_coho$rel.LL <- NULL
out_riv_coho$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_anomaly_night_coho=c(fits_anomaly_night_coho,list(fit))

weights <- akaike.weights(out_riv_coho$AICc)

out_riv_coho$deltaAICc <- weights$deltaAIC
out_riv_coho$rel.LL <- weights$rel.LL
out_riv_coho$weights <- weights$weights


min.AICc <- order(out_riv_coho$AICc)
out_riv_coho.ordered <- out_riv_coho[min.AICc, ]
out_riv_coho.ordered

out_riv_coho.ordered$cumulative_weights <- cumsum(out_riv_coho.ordered$weights)

relative_importance_season <- sum(out_riv_coho$weights[out_riv_coho$season==1])
relative_importance_temp_anomaly <- sum(out_riv_coho$weights[out_riv_coho$temp_anomaly==1])
relative_importance_flow_anomaly <- sum(out_riv_coho$weights[out_riv_coho$flow_anomaly==1])
relative_importance_hatchery <- sum(out_riv_coho$weights[out_riv_coho$hatchery==1])
relative_importance_temp_diff <- sum(out_riv_coho$weights[out_riv_coho$temp_diff==1])
relative_importance_flow_diff <- sum(out_riv_coho$weights[out_riv_coho$flow_diff==1])
riv_anomaly_coho <- data.frame(variable = c("season",
                                            "temperature anomaly",
                                            "flow anomaly",
                                            "temperature difference",
                                            "flow difference",
                                            "hatchery"),
                               relative_importance = c(relative_importance_season,
                                                       relative_importance_temp_anomaly,
                                                       relative_importance_flow_anomaly,
                                                       relative_importance_temp_diff,
                                                       relative_importance_flow_diff,
                                                       relative_importance_hatchery))
#round to 2 decimal places

riv_anomaly_coho$relative_importance <- round(riv_anomaly_coho$relative_importance, 2)

# order in decreasing order

riv_anomaly_coho <- riv_anomaly_coho[order(riv_anomaly_coho$relative_importance, decreasing = TRUE),]

riv_anomaly_coho

#save
write.csv(riv_anomaly_coho, 
          file = here("dungeness",
                      "output","dungenss_coho_relative_importance_anomaly_night.csv"))

#plot



predict_coho1_dungeness_anomaly_night <- predict(dungeness_coho_best_model_anomaly_night, type = "ytT", interval = "confidence")

glimpse(predict_coho1_dungeness_anomaly_night$pred)


head(predict_coho1_dungeness_anomaly_night$pred)

predict_coho1_dungeness_anomaly_night$pred$trap <- 'screw'

predict_coho1_dungeness_anomaly_night$pred$year <-  as.numeric(substr(predict_coho1_dungeness_anomaly_night$pred$.rownames, 
                                                                1, 4))

predict_coho1_dungeness_anomaly_night$pred$daynight_category <- ifelse(substr(predict_coho1_dungeness_anomaly_night$pred$.rownames, 
                                                                        6, nchar(predict_coho1_dungeness_anomaly_night$pred$.rownames)) == 'Day', 'day', 'night')

predict_coho1_dungeness_anomaly_night$pred$doy <- predict_coho1_dungeness_anomaly_night$pred$t+120

dungeness_covariates_coho1 <- as.data.frame(t(covariates_coho1_anomaly_night))

dungeness_covariates_coho1$doy <- as.numeric(rownames(dungeness_covariates_coho1))

dungeness_covariates_coho1_long <-  dungeness_covariates_coho1 %>% 
  dplyr::select(doy, starts_with("coho")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_coho1_dungeness_anomaly_night$pred <- predict_coho1_dungeness_anomaly_night$pred %>%
  left_join(dungeness_covariates_coho1_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_coho1_dungeness_anomaly_night$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"), alpha =0.8)+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2, alpha = 0.6)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_coho1_dungeness_anomaly_night$pred, aes(x = doy, y = coho1_hatchery_perhour_diff, 
                                                             color = "hatchery difference"), alpha =0.5) +
  facet_wrap(~year+daynight_category, ncol = 3, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log (coho salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", "hatchery difference" = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(120, 160), breaks = c(130,150))+
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

ggsave(here("dungeness","output","coho1_dungeness_prediction_w_hatchery_anomaly_night.png"), width = 6, height = 8, units = "in", dpi = 300)  


autoplot(dungeness_coho_best_model_anomaly_night, plot.type = "std.model.resids.ytT") +
  geom_point(aes(color = "wild"), size = 1)+
  labs(x = "Day of year")+
  scale_color_manual(name = "", values = c("wild" = "salmon"))+
  theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        strip.background = element_rect(
          color="white", fill="white"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))  

ggsave(here("dungeness","output","coho1_dungeness_residuals_anomaly_night.png"), width = 8, height = 6, units = "in", dpi = 300)  


# second best model (with hatchery)

dungeness_coho_2best_model_anomaly_night <- fits_anomaly_night_coho[[49]]

#save

save(dungeness_coho_2best_model_anomaly_night, file = here("dungeness",
                                                          "output","dungeness_coho_2best_model_anomaly_night.RData"))
######### start here
#ci
tidy(dungeness_coho_2best_model_anomaly_night)

autoplot(dungeness_coho_2best_model_anomaly_night)

predict_coho1_dungeness_anomaly_night <- predict(dungeness_coho_2best_model_anomaly_night, type = "ytT", interval = "confidence")

glimpse(predict_coho1_dungeness_anomaly_night$pred)

ggsave(here("dungeness","output","coho1_dungeness_prediction_w_hatchery_anomaly_night_2best.png"), width = 6, height = 8, units = "in", dpi = 300)  


ggsave(here("dungeness","output","coho1_dungeness_prediction_w_hatchery_anomaly_night_2best.png"), width = 6, height = 8, units = "in", dpi = 300)  





