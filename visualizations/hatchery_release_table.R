#goal - make a table with number of hatchery releases a year


library(tidyverse)
library(here)


#read data

dungeness_chinook <- read.csv(here("data","data_backup_pied_piper_r",
                                   "dungeness","unmarked_hatchery_chinook_corrected.csv"),header=TRUE)
dungeness_coho <- read.csv(here("data","data_backup_pied_piper_r",
                                "dungeness","unmarked_hatchery_coho_corrected.csv"),header=TRUE)

puyallup_chinook <- read.csv(here("data","data_backup_pied_piper_r",
                                  "puyallup","unmarked_hatchery_chinook_corrected_new.csv"),header=TRUE)

puyallup_coho <- read.csv(here("data","data_backup_pied_piper_r",
                               "puyallup","unmarked_hatchery_coho_corrected_new.csv"),header=TRUE)

skagit_chinook <- read.csv(here("data","data_backup_pied_piper_r",
                                "skagit","unmarked_hatchery_chinook_corrected.csv"),header=TRUE)

skagit_coho <- read.csv(here("data","data_backup_pied_piper_r",
                             "skagit","unmarked_hatchery_coho_corrected.csv"),header=TRUE)

#make the dataframe long by pivot longer

dungeness_chinook_long <- dungeness_chinook %>% 
  mutate(river = "Dungeness", species = "chinook") %>%
  rename(Date = chinook0_dungeness_Date, unmarked = chinook0_dungeness_unmarked,
         wild = chinook0_dungeness_wild_num, hatchery = chinook0_dungeness_hatchery_num,
         unmarked_prop = chinook0_dungeness_unmarked_prop,
         wild_prop = chinook0_dungeness_wild_prop
  )

dungeness_coho_long <- dungeness_coho %>%
  mutate(river = "Dungeness", species = "coho") %>%
  rename(Date = coho1_dungeness_Date, unmarked = coho1_dungeness_unmarked,
         wild = coho1_dungeness_wild_num, hatchery = coho1_dungeness_hatchery_num,
         unmarked_prop = coho1_dungeness_unmarked_prop,
         wild_prop = coho1_dungeness_wild_prop
  )

puyallup_chinook_long <- puyallup_chinook %>%
  mutate(river = "Puyallup", species = "chinook") %>%
  rename(Date = chinook0_puyallup_Date, unmarked = chinook0_puyallup_unmarked,
         wild = chinook0_puyallup_wild_num, hatchery = chinook0_puyallup_hatchery_num,
         unmarked_prop = chinook0_puyallup_unmarked_prop,
         wild_prop = chinook0_puyallup_wild_prop
  )

puyallup_coho_long <- puyallup_coho %>%
  mutate(river = "Puyallup", species = "coho") %>%
  rename(Date = coho1_puyallup_Date, unmarked = coho1_puyallup_unmarked,
         wild = coho1_puyallup_wild_num, hatchery = coho1_puyallup_hatchery_num,
         unmarked_prop = coho1_puyallup_unmarked_prop,
         wild_prop = coho1_puyallup_wild_prop
  )

skagit_chinook_long <- skagit_chinook %>%
  mutate(river = "Skagit", species = "chinook") %>%
  rename(Date = chinook0_skagit_Date, unmarked = chinook0_skagit_unmarked,
         wild = chinook0_skagit_wild_num, hatchery = chinook0_skagit_hatchery_num,
         unmarked_prop = chinook0_skagit_unmarked_prop,
         wild_prop = chinook0_skagit_wild_prop
  )

skagit_coho_long <- skagit_coho %>%
  mutate(river = "Skagit", species = "coho") %>%
  rename(Date = coho1_skagit_Date, unmarked = coho1_skagit_unmarked,
         wild = coho1_skagit_wild_num, hatchery = coho1_skagit_hatchery_num,
         unmarked_prop = coho1_skagit_unmarked_prop,
         wild_prop = coho1_skagit_wild_prop
  )


unmarked_hatchery <- bind_rows(dungeness_chinook_long, dungeness_coho_long,
                               puyallup_chinook_long, puyallup_coho_long,
                               skagit_chinook_long,
                               skagit_coho_long)

unmarked_hatchery <- unmarked_hatchery %>% 
  mutate(year = year(Date))


unmarked_hatchery_agg <- unmarked_hatchery %>% 
  group_by(year, species, river) %>% 
  summarize(releases = n()) %>% 
  group_by(species, river) %>% 
  summarize(min_releases = min(releases), max_releases = max(releases))


#read puyallup rmis data

puyallup_chinook_rmis <- read.csv(here("data","data_backup_pied_piper_r",
                                       "puyallup","chinook_subset_rmis.csv"),header=TRUE)


puyallup_chinook_rmis_agg <- puyallup_chinook_rmis %>% 
  # mutate(Date = as.Date(date, format = "%Y/%m/%d")) %>% 
  mutate(year = year(date)) %>% 
  group_by(year) %>% 
  summarize(total_released = sum(total)) %>% 
  mutate(river = "Puyallup", species = "chinook")

puyallup_coho_rmis <- read.csv(here("data","data_backup_pied_piper_r",
                                    "puyallup","coho_subset_rmis.csv"),header=TRUE)

puyallup_coho_rmis_agg <- puyallup_coho_rmis %>% 
  # mutate(Date = as.Date(date, format = "%Y/%m/%d")) %>% 
  mutate(year = year(date)) %>% 
  group_by(year) %>% 
  summarize(total_released = sum(total)) %>% 
  mutate(river = "Puyallup", species = "coho")


#read skagit rmis data

# read from https://www.rmis.org/reports/CSV9481.txt

skagit_rmis <- read.csv(here("data","data_backup_pied_piper_r",
                             "skagit","rmis.csv"),header=TRUE)


skagit_rmis_agg <- skagit_rmis %>% 
  #conver 8 digit numeric date to date
  mutate(date = as.Date(paste(substr(last_release_date,0,4), substr(last_release_date,5,6),
                              substr(last_release_date,7,8), sep = "/"), 
                        format = "%Y/%m/%d")) %>%
  mutate(year = year(date), age = year - brood_year) %>% 
  group_by(year,species,age) %>%
  mutate(total = tagged_adclipped + tagged_unclipped + untagged_adclipped + 
           untagged_unclipped+untagged_unknown) %>%
  dplyr::select(date, year, species, age, total)


skagit_chinook <- skagit_rmis_agg %>% 
  filter(species == 1, age == 1, year > 2009, year < 2023) %>% 
  group_by(year) %>% 
  summarize(total_released = sum(total)) %>% 
  mutate(river = "Skagit", species = "chinook")

skagit_coho <- skagit_rmis_agg %>%
  filter(species == 2, age == 2, year > 2009, year < 2023) %>% 
  group_by(year) %>% 
  summarize(total_released = sum(total)) %>% 
  mutate(river = "Skagit", species = "coho")



rmis <- bind_rows(puyallup_chinook_rmis_agg, puyallup_coho_rmis_agg,
                  skagit_chinook, skagit_coho)


#read dungness rmis data

dungeness_rmis <- read.csv(here("data","data_backup_pied_piper_r",
                                "dungeness","subset_rmis_new.csv"),header=TRUE)

dungeness_rmis_agg <- dungeness_rmis %>% 
  mutate(date = last_release_date) %>%
  mutate(year = year(date)) %>% 
  dplyr::select(year,total,species,age)



dungeness_chinook_rmis <- dungeness_rmis_agg %>% 
  filter(species == 1, age == 0) %>% 
  group_by(year) %>% 
  summarize(total_released = sum(total)) %>% 
  mutate(river = "Dungeness", species = "chinook")


dungeness_coho_rmis <- dungeness_rmis_agg %>%
  filter(species == 2, age == 1) %>% 
  group_by(year) %>% 
  summarize(total_released = sum(total)) %>% 
  mutate(river = "Dungeness", species = "coho")


#bind all rows

rmis <- bind_rows(dungeness_chinook_rmis, dungeness_coho_rmis, 
                  puyallup_chinook_rmis_agg, puyallup_coho_rmis_agg,
                  skagit_chinook, skagit_coho)

unmarked_hatchery_agg2 <- unmarked_hatchery %>% 
  group_by(year, species, river) %>% 
  summarize(releases = n())

#combine with unmarked hatchery agg by year, species ,river

rmis_unmarked_hatchery <- left_join(rmis, unmarked_hatchery_agg2, by = c("year","species","river"))

#drop year = 2011, river = Skagit


rmis_unmarked_hatchery <- rmis_unmarked_hatchery %>% 
  filter(!(year == 2011 & river == "Skagit"))

#save csv

write_csv(rmis_unmarked_hatchery, here("data","rmis_hatchery_release_data.csv"))



#hatchery name data

#read puyallup_2004-2021.csv

puyallup_hatchery <- read.csv(here("data","data_backup_pied_piper_r",
                                   "puyallup","puyallup_2004-2021_rmis.csv"),header=TRUE)

#group by year and aggregate all the unique values in the hatchery_location_name column

puyallup_hatchery_agg <- puyallup_hatchery %>% 
  group_by(year) %>% 
  summarize(hatchery = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(river = "Puyallup", species = "chinook")

#remove GREENWATER ACCLIMATION PD and PUYALLUP TRIBAL HATCHERY
puyallup_hatchery_agg_new <- puyallup_hatchery %>% 
  filter(hatchery_location_name != "GREENWATER ACCLIMATION PD" &
           hatchery_location_name != "PUYALLUP TRIBAL HATCHERY") %>%
  group_by(year) %>% 
  summarize(hatchery = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(river = "Puyallup", species = "chinook")


#read rmis.csv in skagit


skagit_hatchery <- read.csv(here("data","data_backup_pied_piper_r",
                                 "skagit","rmis.csv"),header=TRUE)


skagit_hatchery_agg_chinook <- skagit_hatchery %>% 
  #conver 8 digit numeric date to date
  mutate(date = as.Date(paste(substr(last_release_date,0,4), substr(last_release_date,5,6),
                              substr(last_release_date,7,8), sep = "/"), 
                        format = "%Y/%m/%d")) %>%
  mutate(year = year(date), age = year - brood_year) %>% 
  group_by(year,species,age) %>%
  mutate(total = tagged_adclipped + tagged_unclipped + untagged_adclipped + 
           untagged_unclipped+untagged_unknown) %>%
  dplyr::select(date, year, species, age, total, hatchery_location_name) %>% 
  group_by(year,species,age) %>% 
  summarize(hatchery = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(river = "Skagit") %>% 
  filter(species == 1, age == 1, year > 2009, year < 2023) %>% 
  dplyr::select(year, hatchery, river) %>% 
  mutate(species = "chinook")


skagit_hatchery_agg_coho <- skagit_hatchery %>% 
  #conver 8 digit numeric date to date
  mutate(date = as.Date(paste(substr(last_release_date,0,4), substr(last_release_date,5,6),
                              substr(last_release_date,7,8), sep = "/"), 
                        format = "%Y/%m/%d")) %>%
  mutate(year = year(date), age = year - brood_year) %>% 
  group_by(year,species,age) %>%
  mutate(total = tagged_adclipped + tagged_unclipped + untagged_adclipped + 
           untagged_unclipped+untagged_unknown) %>%
  dplyr::select(date, year, species, age, total, hatchery_location_name) %>% 
  group_by(year,species,age) %>% 
  summarize(hatchery = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(river = "Skagit") %>% 
  filter(species == 2, age == 2, year > 2009, year < 2023) %>% 
  dplyr::select(year, hatchery, river) %>% 
  mutate(species = "coho")

dungeness_hatchery <- read.csv(here("data","data_backup_pied_piper_r",
                                    "dungeness","subset_rmis_new.csv"),header=TRUE)


dungeness_hatchery_agg_chinook <- dungeness_rmis %>% 
  mutate(date = last_release_date) %>%
  mutate(year = year(date)) %>% 
  dplyr::select(year,total,species,age, hatchery_location_name) %>% 
  group_by(year,species,age) %>% 
  summarize(hatchery = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(river = "Dungeness") %>% 
  filter(species == 1, age == 0) %>% 
  dplyr::select(year, hatchery, river) %>% 
  mutate(species = "chinook")

dungeness_hatchery_agg_coho <- dungeness_rmis %>% 
  mutate(date = last_release_date) %>%
  mutate(year = year(date)) %>% 
  dplyr::select(year,total,species,age, hatchery_location_name) %>% 
  group_by(year,species,age) %>% 
  summarize(hatchery = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(river = "Dungeness") %>% 
  filter(species == 2, age == 1) %>% 
  dplyr::select(year, hatchery, river) %>% 
  mutate(species = "coho")

hatchery_locations <- rbind(dungeness_hatchery_agg_chinook,
                            dungeness_hatchery_agg_coho,
                            puyallup_hatchery_agg,
                            skagit_hatchery_agg_chinook,
                            skagit_hatchery_agg_coho
                            
)

rmis_hatchery_location <- rmis_unmarked_hatchery  %>% 
  left_join(hatchery_locations, by = c("year","species","river")) %>% 
  select(river,species,year, releases, total_released, hatchery)


rmis_hatchery_location_mean <- rmis_hatchery_location %>% 
  group_by(species,river) %>% 
  summarize(average_releases = round(mean(releases, na.rm = TRUE)),
            average_total_released = round(mean(total_released, na.rm = TRUE)),
            hatcheries = paste(unique(hatchery), collapse = ", "))


#save csv
write.csv(rmis_hatchery_location, 
          here("data","rmis_hatchery_release_data_location.csv"))

write.csv(rmis_hatchery_location_mean, 
          here("data","rmis_hatchery_release_data_location_mean.csv"))