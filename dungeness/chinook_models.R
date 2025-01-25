# Goal - To run marss analysis on chinook data with
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

data_day_night <- rbind(data_day,data_night)

#calculate the cumulative sum of chinook0_wild_perhour

data_day_night <- data_day_night %>%
  filter(doy > 120, doy < 250) %>%
  group_by(year, daytime_category) %>% 
  mutate(chinook0_wild_perhour_cumsum = cumsum(ifelse(is.na(chinook0_wild_perhour),
                                                      0,
                                                      chinook0_wild_perhour))) %>% 
  mutate(chinook0_proportion = chinook0_wild_perhour_cumsum/
           sum(chinook0_wild_perhour,na.rm = TRUE))


#calculate the median day of migration for each year (day of year that 0.5 proportion of fish have passed)

# make column for median day of migration - day of year and call it season

median <- data_day_night %>% 
  group_by(year, daytime_category) %>% 
  select(year,chinook0_proportion,doy) %>%
  summarise(median_doy = doy[which.min(abs(chinook0_proportion - 0.5))])

data_day_night <- left_join(data_day_night,median,by = c("year","daytime_category"))

#calculate the difference between the median day of migration and the day of year

data_day_night <- data_day_night %>%
  mutate(season = median_doy - doy)

#plot doy and season, photo_diff, temp_diff, flow_diff, flow, lunar phase, resid

data_day_night %>%
  filter(daytime_category == "Day") %>%
  ggplot() +
  # geom_line(aes(x = doy, y = scale(photoperiod), group = year), col = "darkred") +
  geom_line(aes(x = doy, y = season, group = year), col = "darkblue") +
  geom_line(aes(x = doy, y = scale(season, center = FALSE), group = year), col = "darkred") +
  ylim(-1,1)+
  # geom_line(aes(x = doy, y = photo_diff, group = year), col = "darkblue") +
  # geom_line(aes(x = doy, y = temp_diff, group = year), col = "darkgreen") +
  # geom_line(aes(x = doy, y = flow_diff, group = year), col = "darkorange") +
  # geom_line(aes(x = doy, y = flow, group = year), col = "darkviolet") +
  # geom_line(aes(x = doy, y = lunar_phase, group = year), col = "darkgrey") +
  # geom_line(aes(x = doy, y = resid, group = year), col = "black") +
  theme_minimal()

data_day_night %>% 
  ungroup() %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  select(flow,lunar_phase, season, resid, temp_diff, flow_diff, photo_diff) %>%
  GGally::ggpairs(aes(alpha = 0.2))

ggsave(here("dungeness","output","chinook_covariates_correlation.png"),width = 10, height = 10)

covariates_chinook0 <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy >130 & doy <= 200) %>%
  dplyr::select(year,doy, daytime_category, flow,
                lunar_phase, season, 
                resid, temp_diff, flow_diff, photo_diff, 
                chinook0_hatchery_perhour_interpolate) %>%
  pivot_wider(names_from = c(year, daytime_category), 
              values_from = c(flow, lunar_phase, season, 
                              resid, temp_diff, flow_diff, photo_diff, 
                              chinook0_hatchery_perhour_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

num_years = 2020-2005
num_rows = num_years*2
total_covariates = dim(covariates_chinook0)[1]
covariate_num = total_covariates/num_rows

#scale and center flow
for(i in (1):(num_rows*2)){
  if(sum(covariates_chinook0[i,]) != 0){
    print(rownames(covariates_chinook0)[i])
    covariates_chinook0[i,] = scale(covariates_chinook0[i,],
                                    center = TRUE, scale= TRUE)[,1]
  }
}

for(i in (num_rows*2 + 1):(total_covariates)){
  if(sum(covariates_chinook0[i,]) != 0){
    print(rownames(covariates_chinook0)[i])
    covariates_chinook0[i,] = scale(covariates_chinook0[i,],
    center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_chinook_summer_perhour <- arrange(data_day_night,doy) %>%
  filter(year != 2015 & doy > 130 & doy <= 200) %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category) %>%
  pivot_wider(names_from = c(year, daytime_category), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour)[1]){
  subset_chinook_summer_perhour[i,] = scale(subset_chinook_summer_perhour[i,])[,1]
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
write.csv(df_errors, file = here("dungeness",
                                 "output","error_structure_dungenss_chinook.csv"))

out.tab.all.years.unequalq <- NULL
fits.all.years.unequalq <- NULL
nyears = num_years*2
for(kk in c(1,3,4,7)){
  c = covariates_chinook0[((1+(kk-1)*nyears):(kk*nyears)),]
  name_long = rownames(covariates_chinook0)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-9)
  print(name_individual)
  fit.model = c(list(c= c), mod_list(nyears,1,0,FALSE,TRUE))
  fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  ci = tidy(fit)
  out=data.frame(c=name_individual, Estimate = ci[34,2],
                 conf_low = ci[34,4], conf_high = ci[34,5],
                 AICc=fit$AICc)
  out.tab.all.years.unequalq=rbind(out.tab.all.years.unequalq,out)
  fits.all.years.unequalq=c(fits.all.years.unequalq,list(fit))
}

out.tab.all.years.unequalq$delta_AICc <- out.tab.all.years.unequalq$AICc  - min(out.tab.all.years.unequalq$AICc)

out.tab.all.years.unequalq <- out.tab.all.years.unequalq[order(out.tab.all.years.unequalq$AICc),]

write.csv(out.tab.all.years.unequalq, file = here("dungeness",
                                 "output","correlated_covariates_dungenss_chinook.csv"))

#making combination of covariates

get_covariate_combinations <- function(covariates) {
  n <- length(covariates)
  combinations <- lapply(1:n, function(x) combn(covariates, x, simplify = FALSE))
  unlist(combinations, recursive = FALSE)
}

num_years = 15 
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:6)
out.tab_season<- NULL
out_riv <- NULL
fits_season <- list()

amount_done <- length(fits_season)
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
      k = 2
      lunar_phase = 1
    }
    else if(j==3){
      k = 3
      season = 1
    }
    else if(j==4){
      k = 5
      temperature_difference = 1
    }
    
    else if(j==5){
      k = 6
      flow_difference = 1
    }
    else if(j==6){
      k = 8
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-9))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==8){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, TRUE))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery, FALSE, TRUE))
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
  out_riv <- rbind(out,out_riv)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    temperature_difference = ifelse(temperature_difference==1, fit$par$U[which(covariates == 4),], NA),
                    flow = ifelse(flow == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_difference = ifelse(flow_difference == 1, fit$par$U[which(covariates == 5),], NA),
                    lunar_phase = ifelse(lunar_phase == 1, fit$par$U[which(covariates == 2),], NA),
                    hatchery_day = ifelse(hatchery == 1, fit$par$U["day",], NA),
                    hatchery_night = ifelse(hatchery == 1, fit$par$U["night",], NA),
                    AICc=fit$AICc)
  out.tab_season=rbind(out.tab_season,out2)
  fits_season=c(fits_season,list(fit))
  
  
}

out.tab_season$deltaAICc <- out.tab_season$AICc - min(out.tab_season$AICc)
min.AICc <- order(out.tab_season$AICc)
out.tab_season.ordered <- out.tab_season[min.AICc, ]
out.tab_season.ordered

#rounding the estimates

out.tab_season.ordered<- round(out.tab_season.ordered, 2)
#renaming columns to capitalize
colnames(out.tab_season.ordered) <- c("Season", "Temperature Difference", "Flow", 
                                      "Flow Difference", "Lunar Phase", 
                                      "Hatchery, day", "Hatchery, night",
                                      "AICc", "deltaAICc")
#drop AICc
out.tab_season.ordered <- out.tab_season.ordered[, -8]

write.csv(out.tab_season.ordered, 
          file = here("dungeness",
                      "output","model_selection_dungenss_chinook.csv"))


#relative importance

fit.model = c(mod_list(num_rows,0,0, FALSE, TRUE))
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
out_riv <- rbind(out,out_riv)

# out_riv$deltaAICc <- NULL
out_riv$rel.LL <- NULL
out_riv$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_season=c(fits_season,list(fit))

weights <- akaike.weights(out_riv$AICc)

out_riv$deltaAICc <- weights$deltaAIC
out_riv$rel.LL <- weights$rel.LL
out_riv$weights <- weights$weights


min.AICc <- order(out_riv$AICc)
out_riv.ordered <- out_riv[min.AICc, ]
out_riv.ordered

out_riv.ordered$cumulative_weights <- cumsum(out_riv.ordered$weights)

relative_importance_season <- sum(out_riv$weights[out_riv$season==1])
relative_importance_temperature_difference <- sum(out_riv$weights[out_riv$temperature_difference==1])
relative_importance_flow <- sum(out_riv$weights[out_riv$flow==1])
relative_importance_flow_difference <- sum(out_riv$weights[out_riv$flow_difference==1])
relative_importance_lunar_phase <- sum(out_riv$weights[out_riv$lunar_phase==1])
relative_importance_hatchery <- sum(out_riv$weights[out_riv$hatchery==1])

riv_season <- data.frame(variable = c("season",
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

#round the relative importance

riv_season$relative_importance <- round(riv_season$relative_importance, 2)
#order based on decreasing relative importance
riv_season <- riv_season[order(-riv_season$relative_importance),]

#save
write.csv(riv_season, 
          file = here("dungeness",
                      "output","dungenss_chinook_relative_importance.csv"))

#best model

dungeness_chinook_best_model <- fits_season[[which.min(sapply(fits_season, function(x) x$AICc))]]

#save

save(dungeness_chinook_best_model, file = here("dungeness",
                                               "output","dungeness_chinook_best_model.RData"))

#ci
tidy(dungeness_chinook_best_model)

autoplot(dungeness_chinook_best_model)



# Figures

#####




predict_chinook0_dungeness <- predict(dungeness_chinook_best_model, type = "ytT", interval = "confidence")

glimpse(predict_chinook0_dungeness$pred)

head(predict_chinook0_dungeness$pred)

predict_chinook0_dungeness$pred$trap <- 'screw'

predict_chinook0_dungeness$pred$year <-  as.numeric(substr(predict_chinook0_dungeness$pred$.rownames, 
                                                           1, 4))

predict_chinook0_dungeness$pred$daynight_category <- ifelse(substr(predict_chinook0_dungeness$pred$.rownames, 
                                                                   6, nchar(predict_chinook0_dungeness$pred$.rownames)) == 'Day', 'day', 'night')

predict_chinook0_dungeness$pred$doy <- predict_chinook0_dungeness$pred$t+130

dungeness_covariates_chinook0 <- as.data.frame(t(covariates_chinook0))

dungeness_covariates_chinook0$doy <- as.numeric(rownames(dungeness_covariates_chinook0))

dungeness_covariates_chinook0_long <-  dungeness_covariates_chinook0 %>% 
  select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year","daynight_category"),
               names_pattern = "(.*)_(.{4})_(.*)") %>%
  mutate(year = as.numeric(year), trap = 'screw', daynight_category = ifelse(daynight_category == 'Day', 'day', 'night'))

predict_chinook0_dungeness$pred <- predict_chinook0_dungeness$pred %>%
  left_join(dungeness_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_chinook0_dungeness$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_dungeness$pred, aes(x = doy, y = log(chinook0_hatchery_perhour_interpolate+1), 
                                                        color = "hatchery")) +
  facet_wrap(~year+daynight_category, ncol = 5, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
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

ggsave(here("dungeness","output","chinook0_dungeness_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)  


autoplot(dungeness_chinook_best_model, plot.type = "std.model.resids.ytT") +
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

ggsave(here("dungeness","output","chinook0_dungeness_residuals.png"), width = 8, height = 6, units = "in", dpi = 300)  



data_day_night %>% 
  filter(doy>120, doy <= 200, daytime_category == "Night") %>% 
  ggplot() +
  geom_line(aes(x = doy, y = resid, group = as.factor(year), col = "temperature \n residuals"),alpha = 0.5) +
  geom_line(aes(x = doy, y = temp, group= as.factor(year), col = "temperature"),alpha = 0.5)+
  geom_line(aes(x = doy, y = photoperiod, group = as.factor(year), col = "photoperiod"), 
            linewidth = 2,alpha = 0.7)+
  scale_color_manual(name = "",values = c("temperature \n residuals" = "slategray",
                                          "temperature" = "#B17A79",
                                          "photoperiod" = "#E8D68A")) +
  theme_classic() +
  labs(x = "Day of Year",
       y = "Value",
       title = "") +
  theme(legend.position = "right")

ggsave(here("dungeness","output","dungeness_temperature_residuals.png"), width = 7, height = 5, units = "in", dpi = 300)  


#####
