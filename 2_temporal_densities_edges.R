### Filipe C. Serrano 2025 ###

############################################################################
#                                                                          #
#                  GEOMETRY OF DECLINE: VERTEBRATES                        #
#                                                                          #
############################################################################

library(sp)
library(sf)
library(terra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(sjPlot)
library(DHARMa)
library(MuMIn)
library(plyr)
library(ggstats)
library(performance)
library(lme4)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### IMPORTING DATA ####
trends_populations = read.csv("trends_populations.csv", row.names=NULL) %>% 
  dplyr::filter(Group != "Reptiles") %>% 
  dplyr::mutate(lat_abs = abs(Latitude),
                lat_abs_scaled = scale(lat_abs),
                reldist_trailing_scaled = scale(reldist_trailing),
                reldist_leading_scaled = scale(reldist_leading))



library("rnaturalearth")
world <- map_data("world")

#ideally change to a different crs
ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), colour="black",fill="grey", alpha=.3) + 
  theme_classic(base_size = 14) +
  geom_point(data=trends_populations, aes(x=Longitude, y=Latitude, colour=status), pch=20, size=2.3, alpha = 0.6) +
  scale_color_manual(values = c("black", "darkred", "steelblue4")) +
  labs( x = "Longitude", y = "Latitude", colour = "Population trend") +
  coord_sf(ylim=c(-55,90)) + theme(legend.position = c(.15,.22))

ggsave("map_popstatus.jpeg")


### MODELLING ####

#### DECREASE ####

df_decrease <- trends_populations %>% 
  dplyr::filter(cross_equator == 0)

## 
global_model_decrease <- glmmTMB(formula = decreasing ~ -1 + reldist_trailing_scaled + (reldist_trailing_scaled * Group) + 
                                   (reldist_trailing_scaled * lat_abs_scaled) +
                                   (lat_abs_scaled * Group) +
                                   (reldist_trailing_scaled|Species:Group) +
                                   (1|last_year),
                                 family = binomial,
                                 data = df_decrease,
                                 na.action = "na.fail")


summary(global_model_decrease)
# dredge model selection
dredge_selection_decrease <- dredge(global.model = global_model_decrease)
dredge_selection_decrease
