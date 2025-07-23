# Goal - To run marss analysis on chinook data from the Skagit River with
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


data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In


data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                          na.rm = FALSE)


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)

data_day_night <- data_day_night %>%
  filter(doy > 20, doy < 200) %>%
  group_by(year, daytime_category) %>% 
  mutate(chinook0_wild_perhour_cumsum = cumsum(ifelse(is.na(chinook0_wild_perhour),
                                                      0,
                                                      chinook0_wild_perhour))) %>% 
  mutate(chinook0_proportion = chinook0_wild_perhour_cumsum/
           sum(chinook0_wild_perhour,na.rm = TRUE))

median <- data_day_night %>% 
  group_by(year, daytime_category) %>% 
  select(year,chinook0_proportion,doy) %>%
  summarise(median_doy = doy[which.min(abs(chinook0_proportion - 0.5))])

data_day_night <- left_join(data_day_night,median,by = c("year","daytime_category"))


#calculate the difference between the median day of migration and the day of year

data_day_night <- data_day_night %>%
  mutate(season = median_doy - doy)


data_day_night %>% 
  ungroup() %>%
  filter(doy> 150 & doy <= 189) %>%
  select(flow,lunar_phase, season, resid, temp_diff, flow_diff, photo_diff) %>%
  GGally::ggpairs(aes(alpha = 0.2))

ggsave(here("skagit","output","chinook_covariates_correlation.png"),width = 10, height = 10)



covariates_chinook0_skagit_night <- arrange(data_day_night,doy) %>%
  filter(doy >150 & doy < 189 & daytime_category == 'night') %>%
  dplyr::select(year,doy, daytime_category, flow,
                lunar_phase, season, 
                resid, temp_diff, flow_diff, photo_diff,
                chinook0_hatchery_perhour_inp, trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = c(
    flow,  lunar_phase,season,
    resid, temp_diff, flow_diff, photo_diff, chinook0_hatchery_perhour_inp), names_vary = "fastest") %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years*2
# num_covariates = 10
total_covariates = dim(covariates_chinook0_skagit_night)[1]
covariate_num = total_covariates/num_rows

for(i in 1:(num_rows*2)){ # everything except diffs and hatchery
  print(row.names(covariates_chinook0_skagit_night)[i])
  covariates_chinook0_skagit_night[i,] = scale(covariates_chinook0_skagit_night[i,])[,1]
}

#just scale


for(i in (num_rows*2 + 1):(total_covariates)){
  print(row.names(covariates_chinook0_skagit_night)[i])
  if(sum(covariates_chinook0_skagit_night[i,]) != 0){
    covariates_chinook0_skagit_night[i,] = scale(covariates_chinook0_skagit_night[i,],
                                                 center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_chinook_summer_perhour_night <- arrange(data_day_night,doy) %>%
  filter(doy > 150 & doy < 189 & daytime_category == 'night') %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour_night)[1]){
  subset_chinook_summer_perhour_night[i,] = scale(subset_chinook_summer_perhour_night[i,])[,1]
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

#this part doesn't make sense for skagit

# nyears = num_years*2
# c = NULL
# fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,FALSE))
# fit_equal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
#                      control=list(maxit=2000))
# fit_equal_q$AICc
# 
# fit.model = c(list(c= c), mod_list(nyears,0,0,FALSE,TRUE))
# fit_unequal_q <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
#                        control=list(maxit=2000))
# fit_unequal_q$AICc
# 
# df_errors <- data.frame(Error_structure = c("Equal","Unequal"), 
#                         AICc = c(fit_equal_q$AICc,fit_unequal_q$AICc))

num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:8)
out.tab_season<- NULL
out_riv <- NULL
fits_season <- list()

number_done = length(fits_season)
for(i in (1+number_done):length(list_combinations)){
  
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
  photo_difference = 0
  resid = 0
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
      k = 4
      resid = 1
    }
    else if(j==7){
      k = 7
      photo_difference = 1
    }
    else if(j==8){
      k = 8
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_skagit_night[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_skagit_night)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-17))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  
  if(k==8){
    has_hatchery = 1
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  }
  else{
    has_hatchery = 0
    c_num <- length(covariates)
    fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  }
  fit <- MARSS(subset_chinook_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season,
                 temperature_difference = temperature_difference, flow = flow,
                 flow_difference = flow_difference, resid = resid,
                 photo_difference = photo_difference,
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
                    resid = ifelse(resid == 1, fit$par$U[which(covariates == 6),], NA),
                    photo_difference = ifelse(photo_difference == 1, fit$par$U[which(covariates == 7),], NA),
                    hatchery = ifelse(hatchery == 1, fit$par$U[which(covariates == 8),], NA),
                    AICc=fit$AICc)
  out.tab_season=rbind(out.tab_season,out2)
  fits_season=c(fits_season,list(fit))
  
  
}

out.tab_season$deltaAICc <- out.tab_season$AICc - min(out.tab_season$AICc)
min.AICc <- order(out.tab_season$AICc)
out.tab_season.ordered <- out.tab_season[min.AICc, ]
out.tab_season.ordered


write.csv(out.tab_season.ordered, 
          file = here("skagit",
                      "output","model_selection_skagit_chinook.csv"))



#relative importance

fit.model = c(mod_list(num_rows,0,0))
fit <- MARSS(subset_chinook_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temperature_difference = 0
flow = 0
lunar_phase = 0
hatchery = 0
flow_difference = 0
resid = 0
photo_difference = 0

out=data.frame(c=name, season = season,
               temperature_difference = temperature_difference, flow = flow,
               flow_difference = flow_difference,
               lunar_phase = lunar_phase, 
               resid = resid, photo_difference = photo_difference,
               hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out_riv <- rbind(out,out_riv)

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
relative_importance_resid <- sum(out_riv$weights[out_riv$resid==1])
relative_importance_photo_difference <- sum(out_riv$weights[out_riv$photo_difference==1])

riv_season <- data.frame(variable = c("season",
                                      "temperature difference",
                                      "flow",
                                      "flow difference",
                                      "lunar phase",
                                      "residuals",
                                      "photoperiod difference",
                                      "hatchery"),
                         relative_importance = c(relative_importance_season,
                                                 relative_importance_temperature_difference,
                                                 relative_importance_flow,
                                                 relative_importance_flow_difference,
                                                 relative_importance_lunar_phase,
                                                 relative_importance_resid,
                                                 relative_importance_photo_difference,
                                                 relative_importance_hatchery))


#save
write.csv(riv_season[order(riv_season$relative_importance, decreasing = TRUE),], 
          file = here("skagit",
                      "output","skagit_chinook_relative_importance.csv"))



#best model

skagit_chinook_best_model <- fits_season[[which.min(sapply(fits_season, function(x) x$AICc))]]

#save

save(skagit_chinook_best_model, file = here("skagit",
                                               "output","skagit_chinook_best_model.RData"))

#ci
tidy(skagit_chinook_best_model)

autoplot(skagit_chinook_best_model)


predict_chinook0_skagit <- predict(skagit_chinook_best_model, type = "ytT", interval = "confidence")
glimpse(predict_chinook0_skagit)

head(predict_chinook0_skagit$pred)


predict_chinook0_skagit$pred$trap <- substring(predict_chinook0_skagit$pred$.rownames, 12, 
                                            length(predict_chinook0_skagit$pred$.rownames))

predict_chinook0_skagit$pred$year <- substring(predict_chinook0_skagit$pred$.rownames, 1,4)

predict_chinook0_skagit$pred$daynight_category <- substring(predict_chinook0_skagit$pred$.rownames, 6,10)

skagit_covariates_chinook0 <- as.data.frame(t(covariates_chinook0_skagit_night))
glimpse(skagit_covariates_chinook0)
#make column for doy
skagit_covariates_chinook0$doy <- as.numeric(rownames(skagit_covariates_chinook0))



#make skagit_covariates_chinook0 long by taking year, daynight_category, trap out of the column names

skagit_covariates_chinook0_long <- skagit_covariates_chinook0 %>%
  select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year", "daynight_category","trap"),
               names_pattern = "(.*)_(.{4})_(.{5})_(.{5})")
head(skagit_covariates_chinook0_long)
predict_chinook0_skagit$pred$doy <- predict_chinook0_skagit$pred$t+150

## merge skagit_covariates_chinook0_long$chinook0_hatchery_perhour_inp with predict_chinook0_skagit$pred
# 
predict_chinook0_skagit$pred <- predict_chinook0_skagit$pred %>%
  left_join(skagit_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_chinook0_skagit$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"))+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_skagit$pred, aes(x = doy, y = log(chinook0_hatchery_perhour_inp+1),
                                                    color = "hatchery")) +
  facet_wrap(~year+trap, ncol = 5, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log(Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", hatchery = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(150, 190), breaks = c(130,150,170))+
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

ggsave(here("skagit","output","chinook0_skagit_prediction_w_hatchery.png"), width = 6, height = 8, units = "in", dpi = 300)  



autoplot(skagit_chinook_best_model, plot.type = "std.model.resids.ytT") +
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

ggsave(here("skagit","output","chinook0_skagit_residuals.png"), width = 8, height = 6, units = "in", dpi = 300)  



# anomaly models


######



data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$chinook0_wild_perhour <- data_day_night$chinook0_wild_num/data_day_night$In
data_day_night$chinook0_hatchery_perhour <- data_day_night$chinook0_hatchery_num/data_day_night$In


data_day_night$chinook0_hatchery_perhour_inp <- na.approx(data_day_night$chinook0_hatchery_perhour, 
                                                          na.rm = FALSE)


lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)

data_day_night <- data_day_night %>%
  filter(doy > 20, doy < 200) %>%
  group_by(year, daytime_category) %>% 
  mutate(chinook0_wild_perhour_cumsum = cumsum(ifelse(is.na(chinook0_wild_perhour),
                                                      0,
                                                      chinook0_wild_perhour))) %>% 
  mutate(chinook0_proportion = chinook0_wild_perhour_cumsum/
           sum(chinook0_wild_perhour,na.rm = TRUE)) %>%
  mutate(chinook0_hatchery_perhour_diff = c(NA, diff(chinook0_hatchery_perhour_inp,1)),
         temp_rolling_mean = rollmean(temp, k = 32, fill = NA, align = "right"),
         flow_rolling_mean = rollmean(flow, k = 32, fill = NA, align = "right"),
         temp_anomaly = temp - temp_rolling_mean,
         flow_anomaly = flow - flow_rolling_mean)



data_day_night %>% 
  filter(doy>150, doy <= 200) %>% 
  ggplot() +
  # geom_line(aes(x = doy, y = resid, group = as.factor(year), col = "temperature \n residuals"),alpha = 0.5) +
  # geom_line(aes(x = doy, y = temp, group= as.factor(year), col = "temperature"),alpha = 0.5)+
  # geom_line(aes(x = doy, y = temp_anomaly, group= as.factor(year), col = "temperature \n anomaly"),alpha = 0.5)+
  # geom_line(aes(x = doy, y = photoperiod, group = as.factor(year), col = "photoperiod"),
  # linewidth = 2,alpha = 0.7)+
  # scale_color_manual(name = "",values = c("temperature \n anomaly" = "slategray",
  #                                         "temperature" = "#B17A79")) +
  geom_line(aes(x = doy, y = flow, group= as.factor(year), col = "flow"),alpha = 0.5)+
  geom_line(aes(x = doy, y = flow_anomaly, group= as.factor(year), col = "flow \n anomaly"),alpha = 0.5)+
  scale_color_manual(name = "",values = c("flow \n anomaly" = "slategray",
                                          "flow" = "#93BCCD")) +
  # "photoperiod" = "#E8D68A")) +
  theme_classic() +
  labs(x = "Day of Year",
       y = "Value",
       title = "") +
  theme(legend.position = "right")


ggsave(here("skagit","output","skagit_flow_anomaly.png"), width = 7, height = 5, units = "in", dpi = 300)  



median <- data_day_night %>% 
  group_by(year, daytime_category) %>% 
  dplyr::select(year,chinook0_proportion,doy) %>%
  summarise(median_doy = doy[which.min(abs(chinook0_proportion - 0.5))])

data_day_night <- left_join(data_day_night,median,by = c("year","daytime_category"))


#calculate the difference between the median day of migration and the day of year

data_day_night <- data_day_night %>%
  mutate(season = median_doy - doy)


data_day_night %>% 
  ungroup() %>%
  filter(doy> 150 & doy <= 189) %>%
  dplyr::select(flow_anomaly,
         temp_anomaly,
         # lunar_phase, 
         season, 
         # resid, 
         temp_diff, flow_diff) %>%
  GGally::ggpairs(aes(alpha = 0.2))+
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3)

ggsave(here("skagit","output","chinook_covariates_correlation_anomaly.png"),width = 10, height = 10)



covariates_chinook0_skagit_night <- arrange(data_day_night,doy) %>%
  filter(doy >150 & doy < 189 & daytime_category == 'night') %>%
  dplyr::select(year,doy, daytime_category, 
                flow_anomaly, temp_anomaly,
                # lunar_phase, 
                season, 
                # resid, 
                temp_diff, flow_diff,
                # photo_diff,
                chinook0_hatchery_perhour_diff, trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = c(
    flow_anomaly, temp_anomaly,
    # lunar_phase, 
    season, 
    # resid, 
    temp_diff, flow_diff, 
    # photo_diff,
    chinook0_hatchery_perhour_diff), names_vary = "fastest") %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years*2
# num_covariates = 10
total_covariates = dim(covariates_chinook0_skagit_night)[1]
covariate_num = total_covariates/num_rows

for(i in 1:(num_rows*2)){ # everything except diffs and hatchery
  print(row.names(covariates_chinook0_skagit_night)[i])
  covariates_chinook0_skagit_night[i,] = scale(covariates_chinook0_skagit_night[i,])[,1]
}

#just scale


for(i in (num_rows*2 + 1):(total_covariates)){
  print(row.names(covariates_chinook0_skagit_night)[i])
  if(sum(covariates_chinook0_skagit_night[i,]) != 0){
    covariates_chinook0_skagit_night[i,] = scale(covariates_chinook0_skagit_night[i,],
                                                 center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_chinook_summer_perhour_night <- arrange(data_day_night,doy) %>%
  filter(doy > 150 & doy < 189 & daytime_category == 'night') %>%
  mutate(log.value = log(chinook0_wild_perhour + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

for(i in 1:dim(subset_chinook_summer_perhour_night)[1]){
  subset_chinook_summer_perhour_night[i,] = scale(subset_chinook_summer_perhour_night[i,])[,1]
}

num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:6)
out.tab_anomaly<- NULL
out_riv <- NULL
fits_anomaly <- list()

number_done = length(fits_anomaly)
for(i in (1+number_done):length(list_combinations)){
  
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
  # photo_difference = 0
  # resid = 0
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
      temperature_difference = 1
    }
    
    else if(j==5){
      k = 5
      flow_difference = 1
    }
    # else if(j==6){
    #   k = 4
    #   resid = 1
    # }
    # else if(j==7){
    #   k = 7
    #   photo_difference = 1
    # }
    else if(j==6){
      k = 6
      hatchery = 1
    }
    
    c = rbind(c,covariates_chinook0_skagit_night[((1+(k-1)*num_rows):(k*num_rows)),])
    name_long = rownames(covariates_chinook0_skagit_night)[1+(k-1)*num_rows]
    name = paste(name, substr(name_long,1,nchar(name_long)-17))
    
  }
  # print(c)
  
  print(name)
  c_num <- length(covariates)
  fit.model = c(list(c= c), mod_list(num_rows,c_num,ifelse(k==6,1,0)))
  # if(k==8){
  #   has_hatchery = 1
  #   c_num <- length(covariates)
  #   fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  # }
  # else{
  #   has_hatchery = 0
  #   c_num <- length(covariates)
  #   fit.model = c(list(c= c), mod_list(num_rows,c_num,has_hatchery))
  # }
  fit <- MARSS(subset_chinook_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  
  
  out=data.frame(c=name, season = season,
                 temperature_difference = temperature_difference, 
                 flow_anomaly = flow_anomaly,
                 flow_difference = flow_difference, 
                 temp_anomaly = temp_anomaly,
                 # resid = resid,
                 # photo_difference = photo_difference,
                 # lunar_phase = lunar_phase, 
                 hatchery = hatchery,
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence,
                 stringsAsFactors = FALSE)
  print(out)
  out_riv <- rbind(out,out_riv)
  #save the estimates of the covariates in dataframe
  out2 = data.frame(season = ifelse(season==1, fit$par$U[which(covariates == 3),], NA),
                    temperature_difference = ifelse(temperature_difference==1, fit$par$U[which(covariates == 4),], NA),
                    flow_anomaly = ifelse(flow_anomaly == 1, fit$par$U[which(covariates == 1),], NA),
                    flow_difference = ifelse(flow_difference == 1, fit$par$U[which(covariates == 5),], NA),
                    temp_anomaly = ifelse(temp_anomaly == 1, fit$par$U[which(covariates == 2),], NA),
                    # lunar_phase = ifelse(lunar_phase == 1, fit$par$U[which(covariates == 2),], NA),
                    # resid = ifelse(resid == 1, fit$par$U[which(covariates == 6),], NA),
                    # photo_difference = ifelse(photo_difference == 1, fit$par$U[which(covariates == 7),], NA),
                    hatchery = ifelse(hatchery == 1, fit$par$U[which(covariates == 6),], NA),
                    AICc=fit$AICc)
  out.tab_anomaly=rbind(out.tab_anomaly,out2)
  fits_anomaly=c(fits_anomaly,list(fit))
  
  
}

out.tab_anomaly$deltaAICc <- out.tab_anomaly$AICc - min(out.tab_anomaly$AICc)
min.AICc <- order(out.tab_anomaly$AICc)
out.tab_anomaly.ordered <- out.tab_anomaly[min.AICc, ]

#round to 2 decimal places
out.tab_anomaly.ordered<- round(out.tab_anomaly.ordered,2)
#drop AICC
out.tab_anomaly.ordered <- out.tab_anomaly.ordered %>% 
  dplyr::select(-AICc)
out.tab_anomaly.ordered
write.csv(out.tab_anomaly.ordered, 
          file = here("skagit",
                      "output","model_selection_skagit_chinook_anomaly.csv"))



#relative importance

fit.model = c(mod_list(num_rows,0,0))
fit <- MARSS(subset_chinook_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))
name = "None"
season = 0
temperature_difference = 0
flow_anomaly = 0
temp_anomaly = 0
# lunar_phase = 0
hatchery = 0
flow_difference = 0
# resid = 0
# photo_difference = 0

out=data.frame(c=name, season = season,
               temperature_difference = temperature_difference, 
               flow_anomaly = flow_anomaly,
               flow_difference = flow_difference, 
               temp_anomaly = temp_anomaly,
               # resid = resid,
               # photo_difference = photo_difference,
               # lunar_phase = lunar_phase, 
               hatchery = hatchery,
               logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
               num.iter=fit$numIter, converged=!fit$convergence,
               stringsAsFactors = FALSE)

out_riv <- rbind(out,out_riv)

out_riv$rel.LL <- NULL
out_riv$weights <- NULL

# out_riv=rbind(out_riv,out)
fits_anomaly=c(fits_anomaly,list(fit))

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
relative_importance_flow_anomaly <- sum(out_riv$weights[out_riv$flow_anomaly==1])
relative_importance_temp_anomaly <- sum(out_riv$weights[out_riv$temp_anomaly==1])
relative_importance_flow_difference <- sum(out_riv$weights[out_riv$flow_difference==1])
# relative_importance_lunar_phase <- sum(out_riv$weights[out_riv$lunar_phase==1])
relative_importance_hatchery <- sum(out_riv$weights[out_riv$hatchery==1])
# relative_importance_resid <- sum(out_riv$weights[out_riv$resid==1])
# relative_importance_photo_difference <- sum(out_riv$weights[out_riv$photo_difference==1])

riv_anomaly <- data.frame(variable = c("season",
                                      "temperature difference",
                                      "flow anomaly",
                                      "temp anomaly",
                                      "flow difference",
                                      # "lunar phase",
                                      # "residuals",
                                      # "photoperiod difference",
                                      "hatchery"),
                         relative_importance = c(relative_importance_season,
                                                 relative_importance_temperature_difference,
                                                 relative_importance_flow_anomaly,
                                                 relative_importance_temp_anomaly,
                                                 relative_importance_flow_difference,
                                                 # relative_importance_lunar_phase,
                                                 # relative_importance_resid,
                                                 # relative_importance_photo_difference,
                                                 relative_importance_hatchery))

#round to 2 decimal places
riv_anomaly$relative_importance <- round(riv_anomaly$relative_importance,2)

#save
write.csv(riv_anomaly[order(riv_anomaly$relative_importance, decreasing = TRUE),], 
          file = here("skagit",
                      "output","skagit_chinook_relative_importance_anomaly.csv"))



#best model

skagit_chinook_best_model <- fits_anomaly[[which.min(sapply(fits_anomaly, function(x) x$AICc))]]

#save

save(skagit_chinook_best_model, file = here("skagit",
                                            "output","skagit_chinook_best_model_anomaly.RData"))

#ci
tidy(skagit_chinook_best_model)

autoplot(skagit_chinook_best_model)

# Figures


predict_chinook0_skagit <- predict(skagit_chinook_best_model, type = "ytT", interval = "confidence")
glimpse(predict_chinook0_skagit)

head(predict_chinook0_skagit$pred)


predict_chinook0_skagit$pred$trap <- substring(predict_chinook0_skagit$pred$.rownames, 12, 
                                               length(predict_chinook0_skagit$pred$.rownames))

predict_chinook0_skagit$pred$year <- substring(predict_chinook0_skagit$pred$.rownames, 1,4)

predict_chinook0_skagit$pred$daynight_category <- substring(predict_chinook0_skagit$pred$.rownames, 6,10)

skagit_covariates_chinook0 <- as.data.frame(t(covariates_chinook0_skagit_night))
glimpse(skagit_covariates_chinook0)
#make column for doy
skagit_covariates_chinook0$doy <- as.numeric(rownames(skagit_covariates_chinook0))



#make skagit_covariates_chinook0 long by taking year, daynight_category, trap out of the column names

skagit_covariates_chinook0_long <- skagit_covariates_chinook0 %>%
  dplyr::select(doy, starts_with("chinook")) %>%
  pivot_longer(cols = -c(doy), names_to = c(".value","year", "daynight_category","trap"),
               names_pattern = "(.*)_(.{4})_(.{5})_(.{5})")
head(skagit_covariates_chinook0_long)
predict_chinook0_skagit$pred$doy <- predict_chinook0_skagit$pred$t+150

## merge skagit_covariates_chinook0_long$chinook0_hatchery_perhour_inp with predict_chinook0_skagit$pred
# 
predict_chinook0_skagit$pred <- predict_chinook0_skagit$pred %>%
  left_join(skagit_covariates_chinook0_long, by = c("doy","year", "daynight_category", "trap"))

ggplot(data = predict_chinook0_skagit$pred)+
  
  geom_line(aes(x = doy, y = estimate, color = "wild, predicted"), alpha = 0.8)+
  geom_point(aes(x = doy, y = y, color = "wild, observed"), size = 0.2, alpha = 0.6)+
  geom_ribbon(aes(x = doy, ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.2)+
  geom_line(data = predict_chinook0_skagit$pred, aes(x = doy, y = chinook0_hatchery_perhour_diff,
                                                     color = "hatchery difference", alpha =0.5)) +
  facet_wrap(~year+trap, ncol = 4, labeller = label_wrap_gen(multi_line=FALSE))+
  labs(x = "Day of year", y = "Log (Chinook salmon per hour)", title = "")+
  scale_color_manual(name = "", values = c("wild, predicted" = "salmon", "wild, observed" = "salmon", 
                                           "hatchery difference" = "cadetblue"),
                     guide = guide_legend(override.aes = list(
                       linetype = c(1,NA,1),
                       shape = c(NA,19,NA),
                       size = c(4,2,4))))+
  scale_y_continuous(breaks = c(-3,0,3))+
  scale_x_continuous(limit = c(150, 190), breaks = c(130,150,170))+
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

ggsave(here("skagit","output","chinook0_skagit_prediction_w_hatchery_anomaly.png"), width = 6, height = 8, units = "in", dpi = 300)  



autoplot(skagit_chinook_best_model, plot.type = "std.model.resids.ytT") +
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

ggsave(here("skagit","output","chinook0_skagit_residuals_anomaly.png"), width = 8, height = 6, units = "in", dpi = 300)  


