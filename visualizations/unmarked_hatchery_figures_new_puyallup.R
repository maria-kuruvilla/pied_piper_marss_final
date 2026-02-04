
library(here)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggrepel)
#read files

dungeness_chinook <- read.csv(here("data", "data_backup_pied_piper_r",
                                   "dungeness","unmarked_hatchery_chinook_corrected.csv"),header=TRUE)
dungeness_coho <- read.csv(here("data","data_backup_pied_piper_r",
                                "dungeness","unmarked_hatchery_coho_corrected.csv"),header=TRUE)

puyallup_chinook <- read.csv(here("puyallup", "output",
                                  "unmarked_hatchery_chinook_corrected_jan2026.csv"))

puyallup_coho <- read.csv(here("puyallup", "output",
                               "unmarked_hatchery_coho_corrected_jan2026.csv"))


skagit_chinook <- read.csv(here("data","data_backup_pied_piper_r",
                                "skagit","unmarked_hatchery_chinook_corrected.csv"),header=TRUE)

skagit_coho <- read.csv(here("data","data_backup_pied_piper_r",
                             "skagit","unmarked_hatchery_coho_corrected.csv"),header=TRUE)


#aggregate the data

colnames(dungeness_chinook)
colnames(dungeness_coho)
colnames(puyallup_chinook)
colnames(puyallup_coho)
colnames(skagit_chinook)
colnames(skagit_coho)

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
         # wild_prop = chinook0_puyallup_wild_prop
  )

puyallup_coho_long <- puyallup_coho %>%
  mutate(river = "Puyallup", species = "coho") %>%
  rename(Date = coho1_puyallup_Date, unmarked = coho1_puyallup_unmarked,
         wild = coho1_puyallup_wild_num, hatchery = coho1_puyallup_hatchery_num,
         unmarked_prop = coho1_puyallup_unmarked_prop,
         # wild_prop = coho1_puyallup_wild_prop
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
#caluclate correlation between unmarked and wild for each species and river

correlations <- unmarked_hatchery %>%
  group_by(river, species) %>%
  mutate(max_unmarked = max(unmarked), max_wild = max(wild, na.rm = TRUE)) %>% 
  summarise(correlation = cor(unmarked, wild, use = "complete.obs"), max_unmarked = max(max_unmarked),
            max_wild = max(max_wild))

#also calculate correlations without Puyallup chinook 2019 data and 
#Puyallup coho 2021 data

correlations_new <- unmarked_hatchery %>%
  filter(!(species == "chinook" & river == "Puyallup" & year(Date) == 2019),
         !(species == "coho" & river == "Puyallup" & year(Date) == 2021)
  ) %>%
  group_by(river, species) %>%
  mutate(max_unmarked = max(unmarked), max_wild = max(wild, na.rm = TRUE)) %>% 
  summarise(correlation = cor(unmarked, wild, use = "complete.obs"), max_unmarked = max(max_unmarked),
            max_wild = max(max_wild))


chinook_p <- unmarked_hatchery %>%
  filter(species == "chinook") %>%
  ggplot(aes(x = unmarked, y = wild)) +
  geom_point(alpha = 0.3,size= 2) +
  #label 2019 data in Puyallup
  geom_label_repel(data = unmarked_hatchery  %>% 
              filter(species == "chinook"),# river == "Puyallup",
                     # year(Date) == "2019"),
            aes(label = year(Date)),
            fill = alpha("white", 0.5), # semi-transparent background
            colour = "gray", # solid text color
            label.size = NA, # remove border
            label.padding = 0.15, # padding around text
            size = 3
  ) +
  facet_wrap(~river, ncol = 3, scales = "free") +
  labs(x  = "", y = "Number of\nwild salmon", title = "Chinook") +
  theme_classic() +
  geom_text(data = correlations %>% filter(species == "chinook"), 
            aes(0.03 * max_unmarked,
                0.9 * max_wild,
                label = paste("r =",round(correlation,2)
                )),
            size = 5, color = "black"
  )+
  #include correlation without Puyallup 2019 data
  geom_text(data = correlations_new %>% filter(species == "chinook", river == "Puyallup"),
            aes(3000,
                300,
                label = paste("r (excl. 2019) =",round(correlation,2)
                )),
            size = 5, color = "black"
  )+
  
  #captilize the first letter of the river name
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 16, hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 18)
  ) +
  #make x axis log scale and do not have scientific notation
  scale_x_log10(labels = scales::label_number(accuracy = 1))

chinook_p

coho_p <- unmarked_hatchery %>%
  filter(species == "coho") %>%
  ggplot(aes(x = unmarked, y = wild)) +
  geom_point(alpha = 0.3,size= 2) +
  #label 2019 data in Puyallup
  geom_label_repel(data = unmarked_hatchery %>% 
              filter(species == "coho"), 
                   # river == "Puyallup",
            # year(Date) == "2019"),
            aes(label = year(Date)),
            fill = alpha("white", 0.5), # semi-transparent background
            colour = "gray", # solid text color
            label.size = NA, # remove border
            label.padding = 0.15, size = 3
              
  ) +
  facet_wrap(~river, ncol = 3, scales = "free") +
  labs(x  = "", y = "Number of\nwild salmon", title = "Coho") +
  theme_classic() +
  geom_text(data = correlations %>% filter(species == "coho"), 
            aes(0.05 * max_unmarked,
                0.8 * max_wild,
                label = paste("r =",round(correlation,2)
                )),
            size = 5, color = "black"
  )+
  #include correlation without Puyallup 2019 data
  geom_text(data = correlations_new %>% filter(species == "coho", river == "Puyallup"),
            aes(1000,
                80,
                label = paste("r (excl. 2021) =",round(correlation,2)
                )),
            size = 5, color = "black"
  )+
  
  #captilize the first letter of the river name
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 16, hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        title = element_text(size = 18)
  ) +
  #make x axis log scale and do not have scientific notation
  scale_x_log10(labels = scales::label_number(accuracy = 1))

coho_p

#combine the plots
ggarrange(chinook_p, coho_p, ncol = 1, nrow = 2, 
          # common.legend = TRUE, legend = "bottom",
          heights = c(1,1))

ggsave(here("visualizations","output","unmarked_wild_scatter_new_new_puyallup.png"), width = 10, height = 8, units = "in")

