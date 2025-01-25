#goal - to produce figure that has estimate of effect of covariate in model
# for all three river in each species 

#Chinook dungeness best model

library(here)
library(tidyverse)
library(ggpubr)

load(here("dungeness","output",
          "dungeness_chinook_best_model.RData"))
ci_dungeness_chinook_best_model <- tidy(dungeness_chinook_best_model)
dungeness_chinook_plot <- ggplot(ci_dungeness_chinook_best_model[c(34:38),], 
       aes(x = c("Season", 
                 "Temperature\n difference",
                 "Flow\n difference",
                 "Hatchery,\nday", 
                 "Hatchery,\nnight"),
           y = estimate, 
           ymin = conf.low, 
           ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Dungeness River, Chinook sub yearlings"
  )+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits = c("Temperature\n difference",
               "Hatchery,\nnight","Hatchery,\nday" ,"Season",
               "Flow\n difference"
               
    )) + 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  coord_flip()+
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.31,0.31))+
  geom_rect(aes(xmin = 1.5, xmax = 3.5, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 1, y = 0.2, label = "Dungeness"), size = 10)


dungeness_chinook_plot 


#Puyallup chinook

load(here("puyallup","output",
          "puyallup_chinook_best_model.RData"))

ci_puyallup_chinook_best_model <- tidy(puyallup_chinook_best_model)

puyallup_chinook_plot <- ggplot(ci_puyallup_chinook_best_model[c(39:44),], 
       aes(x = c("Flow",
                 "Lunar Phase",
                 "Season",
                 
                 "Flow\n difference",
                 "Hatchery,\nday", 
                 "Hatchery,\nnight"),
           y = estimate, 
           ymin = conf.low, 
           ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Puyallup River, Chinook sub yearlings"
  )+
  coord_flip()+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits =c("Flow",
              "Lunar Phase",
              "Season",
              "Hatchery,\nnight",
              "Hatchery,\nday", 
              
              "Flow\n difference"))+ 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.31,0.31))+
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 1, y = 0.2, label = "Puyallup"), size = 10)

puyallup_chinook_plot

#Skagit Chinook

load(here("skagit","output",
          "skagit_chinook_best_model.RData"))


ci_skagit_chinook_best_model <- tidy(skagit_chinook_best_model)


skagit_chinook_plot <- ggplot(ci_skagit_chinook_best_model[c(27:29),], 
       aes(x = c("Season",
                 
                 "Temperature\n difference",
                 "Hatchery"),
           y = estimate, 
           ymin = conf.low, 
           ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Skagit River, Chinook sub yearlings"
  )+
  coord_flip()+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits =c("Temperature\n difference",
              "Hatchery",
              "Season"))+ 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.31,0.31))+
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 1, y = 0.2, label = "Skagit"), size = 10)

skagit_chinook_plot

ggpubr::ggarrange(dungeness_chinook_plot, puyallup_chinook_plot, skagit_chinook_plot,
                  labels = c("a", "b", "c"), ncol = 1, nrow = 3, font.label = list(size = 28),
                  common.legend = TRUE, legend = "right",
                  widths = c(1, 1, 1), heights = c(1,1,0.7))

ggsave(here("visualizations",
            "output",
            "chinook_covariates_estimates.jpeg"), width = 16, height = 16)


#Dungeness coho

load(here("dungeness","output",
          "dungeness_coho_best_model_night.RData"))

ci_dungeness_coho_best_model <- tidy(dungeness_coho_best_model_night)


dungeness_coho_plot <- ggplot(ci_dungeness_coho_best_model[c(18:20),], 
       aes(x = c("Flow",
                 "Season",
                 "Flow\n difference"),
           y = estimate, 
           ymin = conf.low, 
           ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Dungeness River, Coho sub yearlings"
  )+
  coord_flip()+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits =c("Flow",
               
               "Flow\n difference",
              "Season"))+ 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.31,0.31))+
  # geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 1, y = 0.2, label = "Dungeness"), size = 10)

dungeness_coho_plot


#Puyallup coho

load(here("puyallup","output",
          "puyallup_coho_best_model.RData"))

ci_puyallup_coho_best_model <- tidy(puyallup_coho_best_model)


puayllup_coho_plot <- ggplot(ci_puyallup_coho_best_model[c(40:45),], 
       aes(x = c("Flow",
                 "Season",
                 "Flow\n difference",
                 "Temperature\n residuals",
                 
                 "Hatchery,\nday",
                 "Hatchery,\nnight"),
           y = estimate, 
           ymin = conf.low, 
           ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Puyallup River, Coho sub yearlings"
  )+
  coord_flip()+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits =c("Flow",
              "Season",
              "Temperature\n residuals",
              "Hatchery,\nnight",
              "Hatchery,\nday",
              
              "Flow\n difference"))+ 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.31,0.31))+
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 1, y = 0.2, label = "Puyallup"), size = 10)

puayllup_coho_plot


#Skagit coho

load(here("skagit","output",
          "skagit_coho_best_model.RData"))


ci_skagit_coho_best_model <- tidy(skagit_coho_best_model)


skagit_coho_plot <- ggplot(ci_skagit_coho_best_model[c(27:31),], 
       aes(x = c("Flow",
                 "Lunar phase",
                 "Temperature\n difference",
                 "Temperature\n residuals",
                 
                 "Hatchery"),
           y = estimate, 
           ymin = conf.low, 
           ymax = conf.up)) +
  geom_pointrange(size = 1, linewidth = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y =""
       # title = "Skagit River, Coho sub yearlings"
  )+
  coord_flip()+
  theme_classic() +
  theme(axis.title.y=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 10, r = 0, b = 0, l = 10)),
        axis.title.x=element_text(size=24, family = "Sans", 
                                  margin = margin(t = 15, r = 0, b = 0, l = 10)),
        axis.text.y = element_text(size = 24, family = "Sans"),
        axis.text.x=element_text(size=24, family = "Sans")) +
  scale_x_discrete(#guide = guide_axis(n.dodge=3),
    limits =c(
      "Temperature\n residuals",
      "Flow",
              "Lunar phase",
              "Temperature\n difference",
              
              "Hatchery"))+ 
  # scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1, 0.2), limits = c(-0.295,0.295))+
  
  scale_y_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), limits = c(-0.31,0.31))+
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -0.295, ymax = 0.295), col = "cadetblue", alpha = 0.0, fill = "cadetblue")+
  geom_text(aes(x = 1, y = 0.2, label = "Skagit"), size = 10)

skagit_coho_plot

ggpubr::ggarrange(dungeness_coho_plot, puayllup_coho_plot, skagit_coho_plot, 
                  labels = c("a", "b", "c"), ncol = 1, nrow = 3, font.label = list(size = 28),
                  common.legend = TRUE, legend = "right",
                  widths = c(1, 1, 1), heights = c(0.7,1,0.9))

ggsave(here("visualizations",
            "output",
            "coho_covariates_estimates.jpeg"), width = 16, height = 16)









