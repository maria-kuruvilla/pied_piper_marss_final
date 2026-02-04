library(here)
library(tidyverse)

#based on Andrew Bergeer's comment that trap only catches fish released from Voights, Cowskull, Rushing water
rmis_puyallup <- read_csv(here("data","puyallup_rmis_og.csv"))
glimpse(rmis_puyallup)

rmis_puyallup_clean <- rmis_puyallup %>% 
  dplyr::select(species, brood_year, last_release_date, release_location_name, 
                hatchery_location_name, release_stage,
                untagged_unclipped, untagged_unknown, tagged_adclipped, 
                tagged_unclipped, untagged_adclipped) %>% 
  #change YYYYMMDD to date
  mutate(date = as.Date(as.character(last_release_date), format = "%Y%m%d"),
         total = untagged_unclipped + untagged_unknown + tagged_unclipped +
           tagged_adclipped + untagged_adclipped,
         release_year = year(date),
         prop_unmarked = (untagged_unclipped + untagged_unknown)/total,
         unmarked = (untagged_unclipped + untagged_unknown)
         ) %>%
  filter(hatchery_location_name != "PUYALLUP TRIBAL HATCHERY", 
         hatchery_location_name != "PUYALLUP HATCHERY", 
         hatchery_location_name != "GREENWATER ACCLIMATION PD",
         release_location_name != "CLARKS CR 10.0027",
         release_location_name != "CLARKS CRK HATCHERY",
         release_year >= 2004, release_year <= 2021) 

rmis_puyallup_chinook <- rmis_puyallup_clean %>% 
  mutate(age = release_year - brood_year) %>% 
  filter(species == 1, age == 1) %>% 
  mutate(release_year = year(date)) %>% 
  group_by(date, release_year) %>% 
  summarize(total_released = sum(total), 
            hatcheries = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  group_by(release_year) %>% 
  summarize(releases = n(),
            average_release = mean(total_released),
            hatcheries = paste(unique(hatcheries), collapse = ", "))
  
rmis_puyallup_coho <- rmis_puyallup_clean %>% 
  mutate(age = release_year - brood_year) %>% 
  filter(species == 2, age == 2)%>% 
  mutate(release_year = year(date)) %>% 
  group_by(date, release_year) %>% 
  summarize(total_released = sum(total), 
            hatcheries = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  group_by(release_year) %>% 
  summarize(releases = n(),
            average_release = mean(total_released),
            hatcheries = paste(unique(hatcheries), collapse = ", "))

rmis_puyallup_summary <- rmis_puyallup_chinook %>% 
  mutate(species = "chinook") %>% 
  rbind(rmis_puyallup_coho %>% mutate(species = "coho")) %>% 
  group_by(species) %>% 
  summarize(total_releases = sum(releases),
            overall_average_number_release = mean(releases),
            overall_average_release = mean(average_release),
            hatcheries = paste(unique(hatcheries), collapse = "; "))

#save

write.csv(rmis_puyallup_summary, here("visualizations", "output","rmis_puyallup_summary.csv") )

# look at max prop

rmis_puyallup_clean %>% 
  group_by(species) %>% 
  summarize(max_prop_unmarked = max(prop_unmarked)) %>% 
  dplyr::select(species, max_prop_unmarked)

#max unmarked prop without 2021 is 7%
rmis_puyallup_clean %>% 
  filter(release_year != 2021, species == 2) %>% 
  # select(release_year, prop_unmarked) %>% 
  # View() %>% 
  summarize(max_prop_unmarked = max(prop_unmarked), 
            mean_prop_unamrked = mean(prop_unmarked)) %>% 
  dplyr::select(max_prop_unmarked, mean_prop_unamrked)


# read regulr puyallup data

puyallup <- read.csv(here("data","puyallup_final.csv")) %>% 
  mutate(Date = as.Date(Date)) %>% 
  left_join(rmis_puyallup_clean %>% 
  mutate(age = release_year - brood_year) %>% 
  filter(species == 2, age == 2)%>% 
  mutate(release_year = year(date)) %>% 
  group_by(date) %>% 
  summarize(total_coho_released = sum(total), 
            prop_unmarked_coho = mean(prop_unmarked),
            coho_hatcheries = paste(unique(hatchery_location_name), collapse = ", ")), join_by(Date == date)) %>% 
  left_join(rmis_puyallup_clean %>% 
              mutate(age = release_year - brood_year) %>% 
              filter(species == 1, age == 1)%>% 
              mutate(release_year = year(date)) %>% 
              group_by(date) %>% 
              summarize(total_chinook_released = sum(total), 
                        prop_unmarked_chinook = mean(prop_unmarked),
                        chinook_hatcheries = paste(unique(hatchery_location_name), collapse = ", ")), join_by(Date == date))

# for each entry in the total coho released, for 10 days after the release, look at the max
# number of hatchery salmon caught



puyallup_coho_long <- rmis_puyallup_clean %>% 
  mutate(age = release_year - brood_year) %>% 
  filter(species == 2, age == 2)%>% 
  mutate(release_year = year(date)) %>% 
  group_by(date) %>% 
  summarize(total_coho_released = sum(total), 
            # prop_unmarked_coho = mean(prop_unmarked),
            unmarked_coho = sum(unmarked),
            coho_hatcheries = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(prop_unmarked_coho = unmarked_coho / total_coho_released)

puyallup_coho_df <- data.frame(date1 = NA, date2 = NA, prop_unmarked_coho = NA,
                               unmarked_coho = NA,
                               max_hatchery,
                               coho_wild_num = NA)
  
for(d in 1:length(puyallup_coho_long$date)){
  # make dataframe
  date1 <- as.Date(puyallup_coho_long$date[d])
  date2 <- as.Date(puyallup_coho_long$date[d]) + 10
  subset_data <- puyallup %>% filter(Date >= date1, Date <= date2)
  max_hatchery <- which.max(subset_data$coho1_hatchery_num)
  # print(subset_data[max_hatchery,])
  if(length(subset_data[max_hatchery,"coho1_wild_num"]) == 0){
    puyallup_coho_df <- puyallup_coho_df 
  } else{
  puyallup_coho_df <- rbind(puyallup_coho_df, data.frame(date1 = as.Date(date1), date2 = as.Date(date2), 
                                     prop_unmarked_coho=puyallup_coho_long$prop_unmarked_coho[d],
                                     unmarked_coho = puyallup_coho_long$unmarked_coho[d],
                                     max_hatchery = subset_data[max_hatchery, "coho1_hatchery_num"],
                                     coho_wild_num  = subset_data[max_hatchery, "coho1_wild_num"]))
  }
}

puyallup_coho_df 
max(puyallup_coho_df$prop_unmarked_coho, na.rm=T)

cor(puyallup_coho_df$unmarked_coho[-33], puyallup_coho_df$coho_wild_num[-33], use = "complete.obs")
#-0.127
cor(puyallup_coho_df$unmarked_coho, puyallup_coho_df$coho_wild_num, use = "complete.obs")
#0.447
ggplot(puyallup_coho_df)+
  geom_point(aes(x = unmarked_coho, y = coho_wild_num))+
  scale_x_log10()





puyallup_chinook_long <- rmis_puyallup_clean %>% 
  mutate(age = release_year - brood_year) %>% 
  filter(species == 1, age == 1)%>% 
  mutate(release_year = year(date)) %>% 
  group_by(date) %>% 
  summarize(total_chinook_released = sum(total), 
            # prop_unmarked_chinook = mean(prop_unmarked),
            unmarked_chinook = sum(unmarked),
            chinook_hatcheries = paste(unique(hatchery_location_name), collapse = ", ")) %>% 
  mutate(prop_unmarked_chinook = unmarked_chinook / total_chinook_released)

puyallup_chinook_df <- data.frame(date1 = NA, date2 = NA, 
                                  prop_unmarked_chinook = NA,
                                  
                                  unmarked_chinook = NA,
                                  max_hatchery,
                                  chinook_wild_num = NA)


for(d in 1:length(puyallup_chinook_long$date)){
  # make dataframe
  date1 <- as.Date(puyallup_chinook_long$date[d])
  date2 <- as.Date(puyallup_chinook_long$date[d]) + 10
  subset_data <- puyallup %>% filter(Date >= date1, Date <= date2)
  max_hatchery <- which.max(subset_data$chinook0_hatchery_num_day + subset_data$chinook0_hatchery_num_night)
  # print(subset_data[max_hatchery,])
  if(length(subset_data[max_hatchery,"chinook0_wild_num_day"]) == 0 || length(subset_data[max_hatchery,"chinook0_wild_num_night"]) == 0){
    puyallup_chinook_df <- puyallup_chinook_df 
  } else{
    puyallup_chinook_df <- rbind(puyallup_chinook_df, data.frame(date1 = as.Date(date1), date2 = as.Date(date2), 
                                                           prop_unmarked_chinook=puyallup_chinook_long$prop_unmarked_chinook[d],
                                                           unmarked_chinook = puyallup_chinook_long$unmarked_chinook[d],
                                                           max_hatchery = subset_data[max_hatchery, "chinook0_hatchery_num_night"] + subset_data[max_hatchery, "chinook0_hatchery_num_day"],
                                                           chinook_wild_num  = subset_data[max_hatchery, "chinook0_wild_num_night"] + subset_data[max_hatchery, "chinook0_wild_num_day"]))
  }
}

puyallup_chinook_df 
max(puyallup_chinook_df$prop_unmarked_chinook, na.rm=T)

cor(puyallup_chinook_df$unmarked_chinook, puyallup_chinook_df$chinook_wild_num, use = "complete.obs")
ggplot(puyallup_chinook_df)+
  geom_point(aes(x = unmarked_chinook, y = chinook_wild_num))+
  geom_text(aes(x = unmarked_chinook, y = chinook_wild_num, label = year(as.Date(puyallup_chinook_df$date1))), hjust=0, vjust=0)+
  scale_x_log10()


cor(puyallup_chinook_df$unmarked_chinook[-24], puyallup_chinook_df$chinook_wild_num[-24], use = "complete.obs")

ggplot(puyallup_chinook_df[-24,])+
  geom_point(aes(x = unmarked_chinook, y = chinook_wild_num))+
  geom_text(aes(x = unmarked_chinook, y = chinook_wild_num, label = year(as.Date(puyallup_chinook_df$date1[-24]))), hjust=0, vjust=0)+
  scale_x_log10()

#save both files

puyallup_chinook_df

puyallup_chinook_df_edited <- puyallup_chinook_df[-1,] %>% 
  rename(chinook0_puyallup_Date = date1,
         chinook0_puyallup_Date_10 = date2,
         chinook0_puyallup_unmarked_prop = prop_unmarked_chinook,
         chinook0_puyallup_unmarked = unmarked_chinook,
         chinook0_puyallup_hatchery_num = max_hatchery,
         chinook0_puyallup_wild_num = chinook_wild_num) %>% 
  mutate(chinook0_puyallup_Date  = as.Date(chinook0_puyallup_Date),
         chinook0_puyallup_Date_10  = as.Date(chinook0_puyallup_Date_10))

write.csv(puyallup_coho_df_edited, here("puyallup", "output","unmarked_hatchery_coho_corrected_jan2026.csv") )

# do same for coho
puyallup_coho_df_edited <- puyallup_coho_df[-1,] %>% 
  rename(coho1_puyallup_Date = date1,
         coho1_puyallup_Date_10 = date2,
         coho1_puyallup_unmarked_prop = prop_unmarked_coho,
         coho1_puyallup_unmarked = unmarked_coho,
         coho1_puyallup_hatchery_num = max_hatchery,
         coho1_puyallup_wild_num = coho_wild_num) %>% 
  mutate(coho1_puyallup_Date  = as.Date(coho1_puyallup_Date),
         coho1_puyallup_Date_10  = as.Date(coho1_puyallup_Date_10))

write.csv(puyallup_chinook_df_edited, here("puyallup", "output","unmarked_hatchery_chinook_corrected_jan2026.csv") )
