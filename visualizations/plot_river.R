library(sf)
library(ggplot2)
library(here)
library(maps)
library(mapdata)

# Read the spatial data
st_layers(here("data","WBD_17_HU2_GDB.gdb"))
st_layers(here("data","NHD_H_Washington_State_GDB.gdb"))
wa <- st_read(here("data","WBD_17_HU2_GDB.gdb"), layer = "WBDHU2")

nhd_hu10 <- st_read(here("data","NHD_H_Washington_State_GDB.gdb"), 
                    layer = "NHDFlowline")


# nhd_hu10_transform <- nhd_hu10 %>% st_transform(crs = 4326)
# 
# nhd_hu10_2 <- read_sf(here("data","NHD_H_Washington_State_GDB.gdb"), 
#                     layer = "NHDFlowline") %>% 
#   st_transform(crs = 4326)

nhd_skagit <- nhd_hu10[grep("Skagit River", nhd_hu10$gnis_name),]
nhd_dungeness <- nhd_hu10[grep("Dungeness River", nhd_hu10$gnis_name),]
nhd_puyallup <- nhd_hu10[grep("Puyallup River", nhd_hu10$gnis_name),]

nhd_gray_wolf <- nhd_hu10[grep("Gray Wolf", nhd_hu10$gnis_name),]
nhd_sauk <- nhd_hu10[grep("Sauk River", nhd_hu10$gnis_name),]
nhd_cascade <- nhd_hu10[grep("Cascade River", nhd_hu10$gnis_name),]
nhd_carbon <- nhd_hu10[grep("Carbon River", nhd_hu10$gnis_name),]
data_huc10 <- st_read(here("data","WBD_17_HU2_GDB.gdb"), layer = "WBDHU10")

data_huc10_dungeness <- data_huc10[data_huc10$huc10 == 1711002003,]
data_huc10_skagit_lower <- data_huc10[data_huc10$huc10 == 1711000702,]
data_huc10_lower_puyallup <- data_huc10[data_huc10$huc10 == 1711001405,]
data_huc10_upper_puyallup <- data_huc10[data_huc10$huc10 == 1711001402,]
data_huc10_carbon_puyallup <- data_huc10[data_huc10$huc10 == 1711001401,]

#from the data_huc10, filter all names that have "Skagit River" in them

data_huc10_skagit <- data_huc10[grep("Skagit River", data_huc10$name),]

data_huc10_sauk <- data_huc10[grep("Sauk River", data_huc10$name),]
data_huc10_cascade <- data_huc10[grep("Cascade River", data_huc10$name),]

data_huc10[grep("Dungeness River", data_huc10$name),]
data_huc10[grep("Grey Wolf", data_huc10$name),]

#HUC12 	171100200307 - Lower Dungeness River
#HUC12 	171100200306 Name 	Middle Dungeness River
#HUC12 	171100200302	Name 	Upper Dungeness River
#HUC12 	171100200301	Name 	Headwaters Dungeness River
#HUC12 	171100200305	Name 	Lower Gray Wolf River
#HUC12 	171100200304	Name 	Upper Gray Wolf River
ggplot() +
  geom_sf(data = wa, color = "grey", fill = NA)+
  geom_sf(data = data_huc10_dungeness, color = "grey", fill = NA)

ggplot() +
  geom_sf(data = data_huc10_dungeness, color = "grey", fill = NA)

state <- map_data("state")
washington <- subset(state, region=="washington")

wa_map <- ggplot() + 
  coord_fixed(1.3) + 
  # geom_polygon(color="black", fill="gray") + 
  geom_polygon(data=washington, mapping=aes(x=long, y=lat, group=group),
               color="black", fill=NA, alpha = 0.7)+
  theme_classic() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+
  geom_sf(data = data_huc10_dungeness, color = "grey", fill = "#6ea599",
          alpha = 0.6)+
  geom_sf(data = data_huc10_skagit, color = "grey", fill = "#8888a2", 
          alpha = 0.6)+
  geom_sf(data = data_huc10_lower_puyallup, color = "grey", fill = "#6e90a5",
          alpha = 0.6)+
  geom_sf(data = data_huc10_upper_puyallup, color = "grey", fill = "#6e90a5",
          alpha = 0.6)+
  geom_sf(data = data_huc10_carbon_puyallup, color = "grey", fill = "#6e90a5",
          alpha = 0.6)+
  geom_sf(data = data_huc10_sauk, color = "grey", fill = "#8888a2",
          alpha = 0.6)+
  geom_sf(data = data_huc10_cascade, color = "grey", fill = "#8888a2",
          alpha = 0.6)+
  geom_sf(data = st_zm(nhd_skagit), color = "#8888a2", alpha = 0.8)+
  geom_sf(data = st_zm(nhd_dungeness), color = "#6ea599", alpha = 0.8)+
  geom_sf(data = st_zm(nhd_puyallup), color = "#6e90a5", alpha = 0.8)+
  geom_sf(data = st_zm(nhd_carbon), color = "#6e90a5", alpha = 0.8)+
  geom_sf(data = st_zm(nhd_gray_wolf), color = "#6ea599", alpha = 0.8)+
  geom_sf(data = st_zm(nhd_sauk), color = "#8888a2", alpha = 0.8)+
  geom_sf(data = st_zm(nhd_cascade), color = "#8888a2", alpha = 0.8)+
  #remove x and y axis lines
  theme(axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #plot points give latitude and longitude 48.445, -122.325
  geom_point(aes(x = -122.325, y = 48.445), color = "black", size = 2, 
             alpha = 0.5)+
  geom_point(aes(x = -123.128, y = 48.140), color = "black", size = 2, 
             alpha = 0.5)+
  geom_point(aes(x = -122.250, y = 47.196), color = "black", size = 2,
             alpha = 0.5)+
  annotate("text",label = "Skagit River", x = -120, y = 48.4, size = 4, 
            color = "#8888a2")+
  annotate("text",label = "Dungeness\n  River", x = -124, y = 48.5, size = 4,
            color = "#6ea599")+
  annotate("text",label = "Puyallup River", x = -121.1, y = 47.4, size = 4,
            color = "#6e90a5")
  
  

wa_map

ggsave(here("output","river_basins_map.png"), 
       wa_map, width = 10, height = 10, dpi = 300)

library(patchwork)
library(jpeg)
library(magick)

#read image
img <- readJPEG(here("..","pied_piper_MARSS","output",
                    "coho_hatchery_wild_cropped.jpg"))

img_magick <- image_read(here("..","pied_piper_MARSS","output",
                     "coho_hatchery_wild_cropped.jpg")) %>% 
  image_ggplot()

img_magick2 <- image_read(here("..","pied_piper_MARSS","output",
                              "coho_hatchery_wild_cropped_edited.jpg")) %>% 
  image_ggplot()
  
wa_map_read <- image_read(here("..","pied_piper_MARSS","output",
                             "river_basins_map_cropped.png")) %>% 
  image_ggplot()

#make img_magick2 half the size as wa_map and put them together


wa_map / img_magick2 + 

ggsave(here("output","river_basins_map_with_coho.png"), 
       wa_map / img_magick, width = 10, height = 14, dpi = 300)

fig1 <- image_read(here("..","pied_piper_r","output",
                     "chinook_coho_plot_new.png")) %>% 
  image_ggplot()

fig1

wa_map_read/img_magick | fig1

ggsave(here("output","river_basins_map_with_coho_fig1.png"), 
       wa_map_read/img_magick | fig1, width = 10, height = 10, dpi = 300, 
       units = "in")

trial <- (wa_map_read / img_magick2) | fig1

ggsave(here("output","river_basins_map_with_coho_fig1.png"), 
       trial, width = 10, height = 14, dpi = 300)


trial2 <- wa_map  +                  # Add plots on top of each other
  inset_element(img_magick2, left = 0.5, bottom = 0, 
                right = 1, top = 0.4)

ggsave(here("output","river_basins_map_with_coho_trial.png"), 
       trial2, width = 10, height = 14, dpi = 300)

trial_together <- trial2|fig1

ggsave(here("output","river_basins_map_with_coho_fig1_trial.png"), 
       trial_together, width = 10, height = 14, dpi = 300)

fig2 <- image_read(here("..","pied_piper_MARSS","output",
                           "dungeness_chinook_wild_hatchery_all_years_linear.png")) %>%
  image_ggplot()

trial_together <- (trial2/fig2)|fig1

trail11 <- (trial2 + plot_spacer()) / fig2 + plot_layout(heights = c(4, -1.1 ,4.5))

ggsave(here("output","trial11.png"), 
       trail11, width = 10, height = 10, dpi = 300)

ggsave(here("output","river_basins_map_with_coho_fig1_fig2.png"), 
       trial_together, width = 10, height = 10, dpi = 300)

# #let's try huc 12
# 
# data_huc12 <- st_read(here("data","WBD_17_HU2_GDB.gdb"), layer = "WBDHU12")
# lower_dungeness <- data_huc12[data_huc12$huc12 == 171100200307,]
# middle_dungeness <- data_huc12[data_huc12$huc12 == 171100200306,]
# cassalery_creek_dungeness <- data_huc12[data_huc12$huc12 == 171100200204,]
# 
# wa_map2 <- ggplot() + 
#   coord_fixed(1.3) + 
#   # geom_polygon(color="black", fill="gray") + 
#   geom_polygon(data=washington, mapping=aes(x=long, y=lat, group=group),
#                color="black", fill=NA)+
#   theme_classic() + 
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(), axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank())+
#   geom_sf(data = lower_dungeness, color = "grey", fill = "#6ea599")+
#   # geom_sf(data = middle_dungeness, color = "grey", fill = "#6ea599", 
#   #         alpha = 0.5)+
#   geom_sf(data = cassalery_creek_dungeness, color = "grey", fill = "#6ea599")
# 
# wa_map2
 