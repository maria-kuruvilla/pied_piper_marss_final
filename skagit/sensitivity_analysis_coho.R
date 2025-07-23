# goal - to test the sensitivity of the results of Skagit coho population to 
# subtracting x% of hatchery coho from the wild coho numbers

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




data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))

data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In


data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                       na.rm = FALSE)



lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)

data_day_night <- data_day_night %>% add_residuals(lm_temp_day)

percent <- 0.1

data_day_night <- data_day_night %>%
  filter(doy > 20, doy < 200) %>%
  group_by(year, daytime_category) %>% 
  mutate(coho1_wild_perhour_cumsum = cumsum(ifelse(is.na(coho1_wild_perhour),
                                                   0,
                                                   coho1_wild_perhour))) %>% 
  mutate(coho1_proportion = coho1_wild_perhour_cumsum/
           sum(coho1_wild_perhour,na.rm = TRUE)) %>% 
  mutate(coho1_hatchery_perhour_diff = c(NA, diff(coho1_hatchery_perhour_inp,1)),
         temp_rolling_mean = zoo::rollapply(temp, width = 32, FUN = mean, na.rm = TRUE, fill = NA, align = "right"),
         flow_rolling_mean = zoo::rollapply(flow, width = 32,  FUN = mean, na.rm = TRUE, fill = NA, align = "right"),
         temp_anomaly = temp - temp_rolling_mean,
         flow_anomaly = flow - flow_rolling_mean) %>% 
  mutate(percent_wild = ifelse(coho1_wild_perhour - percent*coho1_hatchery_perhour_inp < 0,
                               0,
                               coho1_wild_perhour - percent*coho1_hatchery_perhour_inp),
         percent_hatchery = coho1_hatchery_perhour_inp + percent*coho1_hatchery_perhour_inp) %>%
  mutate(percent_hatchery_diff = c(NA, diff(percent_hatchery,1)),
         coho1_percent_wild_perhour_cumsum = cumsum(ifelse(is.na(percent_wild),
                                                   0,
                                                   percent_wild))) %>%
  mutate(coho1_percent_wild_proportion = coho1_percent_wild_perhour_cumsum/
           sum(percent_wild,na.rm = TRUE))
  
glimpse(data_day_night)

#plot hatchery_diff and percent hatchery_diff ,facet wrap by year

ggplot(data_day_night %>% filter(doy >100, doy < 160, daytime_category == "night"),
       aes(x = doy, y = percent_hatchery_diff)) +
  geom_line() +
  geom_line(aes(y = coho1_hatchery_perhour_diff), color = "red", alpha = 0.5, size=2) +
  facet_wrap(~year) +
  theme_bw() +
  labs(x = "Day of Year", y = "Hatchery Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

median <- data_day_night %>% 
  group_by(year, daytime_category) %>% 
  dplyr::select(year,coho1_percent_wild_proportion,doy) %>%
  summarise(median_doy = doy[which.min(abs(coho1_percent_wild_proportion - 0.5))])

data_day_night <- left_join(data_day_night,median,by = c("year","daytime_category"))



#calculate the difference between the median day of migration and the day of year

data_day_night <- data_day_night %>%
  mutate(season = median_doy - doy)


covariates_coho1_skagit_night <- arrange(data_day_night,doy) %>%
  filter(doy >100 & doy < 150 & daytime_category == 'night') %>%
  dplyr::select(year,doy, daytime_category, flow_anomaly,
                temp_anomaly,
                # lunar_phase, 
                season, 
                # resid, 
                temp_diff, flow_diff, 
                # photo_diff, 
                percent_hatchery_diff, trap) %>% 
  pivot_wider(names_from = c(year, daytime_category, trap), 
              values_from = c(flow_anomaly,
                              temp_anomaly,
                              # lunar_phase, 
                              season, # resid, 
                              temp_diff, flow_diff,  # photo_diff, 
                              percent_hatchery_diff), 
              names_vary = "fastest") %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


#scaling the variables

num_years = 2022-2010
num_rows = num_years*2
# num_covariates = 10
total_covariates = dim(covariates_coho1_skagit_night)[1]
covariate_num = total_covariates/num_rows

for(i in 1:(num_rows*2)){ # everything except diffs and hatchery
  print(row.names(covariates_coho1_skagit_night)[i])
  covariates_coho1_skagit_night[i,] = scale(covariates_coho1_skagit_night[i,])[,1]
}

#just scale


for(i in (num_rows*2 + 1):(total_covariates)){
  print(row.names(covariates_coho1_skagit_night)[i])
  if(sum(covariates_coho1_skagit_night[i,]) != 0){
    covariates_coho1_skagit_night[i,] = scale(covariates_coho1_skagit_night[i,],
                                              center = FALSE, scale= TRUE)[,1]
  }
}

#subset response variable
subset_coho_summer_perhour_night <- arrange(data_day_night,doy) %>%
  filter(doy > 100 & doy < 150 & daytime_category == 'night') %>%
  mutate(log.value = log(percent_wild + 1)) %>%
  dplyr::select(log.value,year,doy,daytime_category,trap) %>%
  pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
  subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
}


#function for C matrix
#######


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



nyears = num_years*2
c<-NULL
for(kk in c(1,3,6)){
  c = rbind(c,covariates_coho1_skagit_night[((1+(kk-1)*nyears):(kk*nyears)),])
  name_long = rownames(covariates_coho1_skagit_night)[1+(kk-1)*nyears]
  name_individual = substr(name_long,1,nchar(name_long)-15)
  print(name_individual)
  
  # out=data.frame(c=name_individual, d = "None",
  #                logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
  #                num.iter=fit$numIter, converged=!fit$convergence,
  #                stringsAsFactors = FALSE)
  # out.tab.all.years_night=rbind(out.tab.all.years_night,out)
  # fits.all.years_night=c(fits.all.years_night,list(fit))
}
fit.model = c(list(c= c), mod_list(nyears,3,1))
fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
             control=list(maxit=2000))

ci_0.1 <- tidy(fit)


#loop through different percent values
percent <- c(0,0.05,0.1,0.15)

sensitivity <- data.frame(matrix(NA, nrow = 3*length(percent), ncol = 3))

sensitivity$variable <- rep(c("flow anomaly","season", "hatchery difference"),length(percent))

sensitivity$percent <- rep(percent, each = 3)





for(p in 1:length(percent)){
  data_day_night <- read.csv(here("skagit", "data", "skagit_2010-2022_w_covariates.csv"))
  
  data_day_night$coho1_wild_perhour <- data_day_night$coho1_wild_num/data_day_night$In
  data_day_night$coho1_hatchery_perhour <- data_day_night$coho1_hatchery_num/data_day_night$In
  
  
  data_day_night$coho1_hatchery_perhour_inp <- na.approx(data_day_night$coho1_hatchery_perhour, 
                                                         na.rm = FALSE)
  
  
  
  lm_temp_day <- lm(temp ~ 1+photoperiod, data = data_day_night)
  
  data_day_night <- data_day_night %>% add_residuals(lm_temp_day)
  
  
  data_day_night <- data_day_night %>%
    filter(doy > 20, doy < 200) %>%
    group_by(year, daytime_category) %>% 
    mutate(coho1_wild_perhour_cumsum = cumsum(ifelse(is.na(coho1_wild_perhour),
                                                     0,
                                                     coho1_wild_perhour))) %>% 
    mutate(coho1_proportion = coho1_wild_perhour_cumsum/
             sum(coho1_wild_perhour,na.rm = TRUE)) %>% 
    mutate(coho1_hatchery_perhour_diff = c(NA, diff(coho1_hatchery_perhour_inp,1)),
           temp_rolling_mean = zoo::rollapply(temp, width = 32, FUN = mean, na.rm = TRUE, fill = NA, align = "right"),
           flow_rolling_mean = zoo::rollapply(flow, width = 32,  FUN = mean, na.rm = TRUE, fill = NA, align = "right"),
           temp_anomaly = temp - temp_rolling_mean,
           flow_anomaly = flow - flow_rolling_mean) %>% 
    mutate(percent_wild = coho1_wild_perhour + percent[p]*coho1_hatchery_perhour_inp,
           percent_hatchery = coho1_hatchery_perhour_inp - percent[p]*coho1_hatchery_perhour_inp) %>%
    mutate(percent_hatchery_diff = c(NA, diff(percent_hatchery,1)),
           coho1_percent_wild_perhour_cumsum = cumsum(ifelse(is.na(percent_wild),
                                                             0,
                                                             percent_wild))) %>%
    mutate(coho1_percent_wild_proportion = coho1_percent_wild_perhour_cumsum/
             sum(percent_wild,na.rm = TRUE))
  
  
  
  
  median <- data_day_night %>% 
    group_by(year, daytime_category) %>% 
    dplyr::select(year,coho1_percent_wild_proportion,doy) %>%
    summarise(median_doy = doy[which.min(abs(coho1_percent_wild_proportion - 0.5))])
  
  data_day_night <- left_join(data_day_night,median,by = c("year","daytime_category"))
  
  
  
  #calculate the difference between the median day of migration and the day of year
  
  data_day_night <- data_day_night %>%
    mutate(season = median_doy - doy)
  
  
  covariates_coho1_skagit_night <- arrange(data_day_night,doy) %>%
    filter(doy >100 & doy < 150 & daytime_category == 'night') %>%
    dplyr::select(year,doy, daytime_category, flow_anomaly,
                  temp_anomaly,
                  # lunar_phase, 
                  season, 
                  # resid, 
                  temp_diff, flow_diff, 
                  # photo_diff, 
                  percent_hatchery_diff, trap) %>% 
    pivot_wider(names_from = c(year, daytime_category, trap), 
                values_from = c(flow_anomaly,
                                temp_anomaly,
                                # lunar_phase, 
                                season, # resid, 
                                temp_diff, flow_diff,  # photo_diff, 
                                percent_hatchery_diff), 
                names_vary = "fastest") %>%
    column_to_rownames(var = "doy") %>%
    as.matrix() %>%
    t()
  
  
  #scaling the variables
  
  num_years = 2022-2010
  num_rows = num_years*2
  # num_covariates = 10
  total_covariates = dim(covariates_coho1_skagit_night)[1]
  covariate_num = total_covariates/num_rows
  
  for(i in 1:(num_rows*2)){ # everything except diffs and hatchery
    print(row.names(covariates_coho1_skagit_night)[i])
    covariates_coho1_skagit_night[i,] = scale(covariates_coho1_skagit_night[i,])[,1]
  }
  
  #just scale
  
  
  for(i in (num_rows*2 + 1):(total_covariates)){
    print(row.names(covariates_coho1_skagit_night)[i])
    if(sum(covariates_coho1_skagit_night[i,]) != 0){
      covariates_coho1_skagit_night[i,] = scale(covariates_coho1_skagit_night[i,],
                                                center = FALSE, scale= TRUE)[,1]
    }
  }
  
  #subset response variable
  subset_coho_summer_perhour_night <- arrange(data_day_night,doy) %>%
    filter(doy > 100 & doy < 150 & daytime_category == 'night') %>%
    mutate(log.value = log(percent_wild + 1)) %>%
    dplyr::select(log.value,year,doy,daytime_category,trap) %>%
    pivot_wider(names_from = c(year, daytime_category, trap), values_from = log.value) %>%
    column_to_rownames(var = "doy") %>%
    as.matrix() %>%
    t()
  
  
  for(i in 1:dim(subset_coho_summer_perhour_night)[1]){
    subset_coho_summer_perhour_night[i,] = scale(subset_coho_summer_perhour_night[i,])[,1]
  }
  nyears = num_years*2
  c<-NULL
  for(kk in c(1,3,6)){
    c = rbind(c,covariates_coho1_skagit_night[((1+(kk-1)*nyears):(kk*nyears)),])
    name_long = rownames(covariates_coho1_skagit_night)[1+(kk-1)*nyears]
    name_individual = substr(name_long,1,nchar(name_long)-15)
    print(name_individual)
    
    # out=data.frame(c=name_individual, d = "None",
    #                logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
    #                num.iter=fit$numIter, converged=!fit$convergence,
    #                stringsAsFactors = FALSE)
    # out.tab.all.years_night=rbind(out.tab.all.years_night,out)
    # fits.all.years_night=c(fits.all.years_night,list(fit))
  }
  fit.model = c(list(c= c), mod_list(nyears,3,1))
  fit <- MARSS(subset_coho_summer_perhour_night, model=fit.model, silent = TRUE, method = "BFGS",
               control=list(maxit=2000))
  ci <- tidy(fit)
  
  sensitivity[(1+3*(p-1)):(3*p),1] <- ci$estimate[27:29]
  sensitivity[(1+3*(p-1)):(3*p),2] <- ci$conf.low[27:29]
  sensitivity[(1+3*(p-1)):(3*p),3] <- ci$conf.up[27:29]
  
}

#rename X1 as estimate
colnames(sensitivity)[1:4] <- c("estimate", "conf.low", "conf.up", "variable")

p4 <- sensitivity %>% 
  filter(variable == "hatchery difference") %>%
  ggplot(aes(x = percent, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up),color = "#725e92") +
  geom_line(color = "#725e92",alpha = 0.8) +
  xlim(-0.2,0.2) +
  # ylim(-0.1,0.1) +
  # geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_smooth(method="lm",fullrange=TRUE,color = "#725e92", alpha = 0.5) +
  labs(#title = "Sensitivity Analysis - Skagit River Coho Salmon",
    x = "Proportion of hatchery coho relabelled as wild coho",
    y = "Estimate of effect of hatchery difference") +
  theme_classic()

p4

ggsave(here("skagit", "output", "sensitivity_analysis_coho.png"),
       plot = p4,
       width = 6,
       height = 4,
       dpi = 300)
