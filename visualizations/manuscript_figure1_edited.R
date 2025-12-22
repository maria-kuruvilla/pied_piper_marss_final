library(sf)
library(ggplot2)
library(here)
library(maps)
library(mapdata)

library(patchwork)
library(jpeg)
library(magick)


state <- map_data("state")
washington <- subset(state, region=="washington")


#rivers
nhd_hu10 <- st_read(here("data","NHD_H_Washington_State_GDB.gdb"), 
                    layer = "NHDFlowline")

nhd_skagit <- nhd_hu10[grep("Skagit River", nhd_hu10$gnis_name),]
nhd_dungeness <- nhd_hu10[grep("Dungeness River", nhd_hu10$gnis_name),]
nhd_puyallup <- nhd_hu10[grep("Puyallup River", nhd_hu10$gnis_name),]
nhd_gray_wolf <- nhd_hu10[grep("Gray Wolf", nhd_hu10$gnis_name),]
nhd_sauk <- nhd_hu10[grep("Sauk River", nhd_hu10$gnis_name),]
nhd_cascade <- nhd_hu10[grep("Cascade River", nhd_hu10$gnis_name),]
nhd_carbon <- nhd_hu10[grep("Carbon River", nhd_hu10$gnis_name),]

#basins
# data_huc10 <- st_read(here("data","WBD_17_HU2_GDB.gdb"), layer = "WBDHU10")

data_huc10 <- st_read(here("..","..","..","Downloads","WBD_17_HU2_GDB",
                           "WBD_17_HU2_GDB.gdb"), 
                      layer = "WBDHU10")

data_huc10_dungeness <- data_huc10[data_huc10$huc10 == 1711002003,]
data_huc10_skagit_lower <- data_huc10[data_huc10$huc10 == 1711000702,]
data_huc10_lower_puyallup <- data_huc10[data_huc10$huc10 == 1711001405,]
data_huc10_upper_puyallup <- data_huc10[data_huc10$huc10 == 1711001402,]
data_huc10_carbon_puyallup <- data_huc10[data_huc10$huc10 == 1711001401,]
data_huc10_skagit <- data_huc10[grep("Skagit River", data_huc10$name),]
data_huc10_sauk <- data_huc10[grep("Sauk River", data_huc10$name),]
data_huc10_cascade <- data_huc10[grep("Cascade River", data_huc10$name),]


wa_map <- ggplot() + 
  coord_fixed(1.3) + 
  # geom_polygon(color="black", fill="gray") + 
  geom_polygon(data=washington, mapping=aes(x=long, y=lat, group=group),
               color="#A1A6AA", fill=NA, 
               linewidth = 0.8)+
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
  #plot points give latitude and longitude 48.445, -122.325
  geom_point(aes(x = -122.325, y = 48.445), color = "black", size = 3, fill = "grey",
             alpha = 0.8, shape = 23)+
  geom_point(aes(x = -123.128, y = 48.140), color = "black", size = 3, fill = "grey",
             alpha = 0.8, shape = 23)+
  geom_point(aes(x = -122.250, y = 47.196), color = "black", size = 3, fill = "grey",
             alpha = 0.8, shape = 23, label = "trap")+
  geom_point(aes(x =-121.735300, y = 48.533900), color = "orange", size = 2,
             alpha = 0.8, shape = 16, label = "coho release site")+
  # geom_point(aes(x = -121.746100, y = 48.433400), color = "#9E6767", size = 2,
             # alpha = 0.8, shape = 15)+
  #plot points of hatchery release locations
  #+48.524200, -121.429200
  geom_point(aes(x = -121.429200, y = 48.524200), color = "#464573", size = 2, 
             alpha = 0.8, shape = 16,label = "chinook release site")+
  #marblemount hatchery +48.433400, -121.746100
  # geom_point(aes(x = -121.746100, y = 48.433400), color = "#9E6767", size = 2, 
  #            alpha = 0.8, shape = 15)+
  #release site for chinook +48.387400, -122.366100
  # geom_point(aes(x = -122.366100, y = 48.387400), color = "#464573", size = 2, 
  #            alpha = 0.8, shape = 16)+
  #chinook release site +48.638400, -121.300000
  geom_point(aes(x = -121.300000, y = 48.638400), color = "#464573", size = 2, 
             alpha = 0.8, shape = 16)+
  #coho release site +48.562400, -121.734100
  # geom_point(aes(x = -121.734100, y = 48.562400), color = "orange", size = 2, 
  #            alpha = 0.8, shape = 16)+
  #chinook release site (+47.976700, -123.110500)
  geom_point(aes(x = -123.110500, y = 47.976700), color = "#464573", size = 2, 
             alpha = 0.8, shape = 16)+
  #coho release site +48.150800, -123.133100
  #cohorelease site 48.028385, -123.139695
  geom_point(aes(x = -123.139695, y = 48.028385), color = "orange", size = 2, 
             alpha = 0.8, shape = 16)+
  #chinook release site +47.199700, -122.257300
  # geom_point(aes(x = -122.257300, y = 47.199700), color = "#464573", size = 2, 
  #            alpha = 0.8, shape = 16)+
  #coho relese site +47.214000, -122.340000
  # geom_point(aes(x = -122.340000, y = 47.214000), color = "orange", size = 2, 
  #            alpha = 0.8, shape = 16)+
  #chinook release site +47.135600, -122.074800
  geom_point(aes(x = -122.074800, y = 47.135600), color = "#464573", size = 2,
             alpha = 0.8, shape = 16)+
  #chinook release site +47.087200, -122.184300
  geom_point(aes(x = -122.184300, y = 47.087200), color = "#464573", size = 2, 
             alpha = 0.8, shape = 16)+
  #coho release site +47.087200, -122.184300
  geom_point(aes(x = -122.184300, y = 47.087200), color = "orange", size = 2, 
             alpha = 0.8, shape = 16)+
  #label
  # geom_point(aes(x = -120.2, y = 47.76), color = "orange", size = 2, 
  #            alpha = 0.8, shape = 16)+
  # #label
  # geom_point(aes(x = -120.2, y = 47.82), color = "#464573", size = 2, 
  #            alpha = 0.8, shape = 16)+
  # 
  # annotate("text",label = "coho", x = -120, y = 47.76, size = 4, 
  #          color = "orange")+
  # annotate("text",label = "chinook", x = -119.9, y = 47.82, 
  #          size = 4, 
  #          color = "#464573")+
  # annotate("text",label = "release sites", x = -120, y = 47.9, 
  #          size = 4, 
  #          color = "slategray")+
  annotate("text",label = "Skagit River", x = -120, y = 48.4, size = 4, 
           color = "#8888a2")+
  annotate("text",label = "Dungeness\n  River", x = -124, y = 48.5, size = 4,
           color = "#6ea599")+
  annotate("text",label = "Puyallup River", x = -121.5, y = 47.4, size = 4,
           color = "#6e90a5")+
  theme_classic()+
  #remove x and y axis lines
  theme(axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
        )


wa_map

# add locations of environmental data collection sites

# 48.1434051232356, -123.128730880501
# Decimal latitude	48.6717922878038	, Decimal longitude	-121.246234788845
# Decimal latitude	47.1851011479701, Decimal longitude	-122.22956251685	
# Latitude: 46 deg; 56 min N
# Longitude:
#   121 deg; 57 min W


wa_map +
  geom_point(aes(x = -123.128730880501, y = 48.1434051232356), 
             color = "darkred", size = 3, shape = 8)+
  geom_point(aes(x = -121.246234788845, y = 48.6717922878038),
             color = "darkred", size = 3, shape = 8)+
  geom_point(aes(x = -122.22956251685, y = 47.1851011479701),
             color = "darkred", size = 3, shape = 8)+
  
  





ggsave(here("output","river_basins_release_map_new_wo_label.png"), 
       wa_map, width = 10, height = 10, dpi = 300)


img_magick2 <- image_read(here("..","pied_piper_MARSS","output",
                               "coho_hatchery_wild_cropped_edited.jpg")) %>% 
  image_ggplot()


trial2 <- wa_map  +                  # Add plots on top of each other
  inset_element(img_magick2, left = 0.55, bottom = 0, 
                right = 1, top = 0.45)+
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm"))
trial2

ggsave(here("output","manuscript_fig1a_release_sites_new.png"), 
       trial2, width = 8, height = 8, units = "in",
       dpi = 300)


fig2 <- image_read(here("..","pied_piper_MARSS","output",
                        "dungeness_chinook_wild_hatchery_all_years_linear_w_hatchery_release.png")) %>%
  image_ggplot()

#make the annotation tag level bold "a", "b", "c"
#remove space above the plot
manuscript_fig1 <- trial2/fig2 + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm"))

ggsave(here("output","manuscript_fig1_release_sites_new.png"), 
       manuscript_fig1, width = 8, height = 8, units = "in",
       dpi = 300)



wa_map2 <- image_read(here("visualizations",
                           "output","river_basins_release_map_new_wo_label_cropped.png")) %>% 
  image_ggplot()

trial3 <- wa_map2  +                  # Add plots on top of each other
  inset_element(img_magick2, left = 0.55, bottom = 0, 
                right = 1, top = 0.37)+
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm"))

manuscript_fig2 <- trial3/fig2 + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm"))


ggsave(here("output","manuscript_fig1_release_sites_crop_new_w_hatchery_release.png"), 
       manuscript_fig2, width = 8, height = 8, units = "in",
       dpi = 300)

#editing the hatchery release

fig3 <- image_read(here("..","pied_piper_MARSS","output",
                        "dungeness_chinook_wild_hatchery_all_years_linear_w_hatchery_release2.png")) %>%
  image_ggplot()


manuscript_fig3 <- trial3/fig3 + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm"))



ggsave(here("output","manuscript_fig1_release_sites_crop_new_w_hatchery_release2.png"), 
       manuscript_fig3, width = 8, height = 8, units = "in",
       dpi = 300)
