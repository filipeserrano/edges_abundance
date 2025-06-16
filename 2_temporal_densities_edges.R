### Filipe C. Serrano 2025 ###

############################################################################
#                                                                          #
#                  GEOMETRY OF DECLINE: VERTEBRATES                        #
#                                                                          #
############################################################################

library(sp)
library(sf)
library(units)
library(terra)a
library(tidyverse)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(mapview)
library(reshape2)
library(lwgeom)
library(ape)
library(adephylo)
library(geiger)
library(readxl)
library(enmSdmX)
library("geosphere")
library(spThin)
library(readr)
library( polylabelr )
library(okara)
library(rasterSp)
library(spatialEco)
library(RColorBrewer)
library(sjPlot)
library(DHARMa)
library(MuMIn)
library(plyr)
library(ggstats)
library(sf)
library(dplyr)
library(units)
library(rlpi)
library("nngeo")
library(rnaturalearthhires)
library(Hmisc)
library(stars)
library(performance)
library(glmmTMB)
library(ggtext)
library(ggeffects)

"%nin% " <- function(x, table) match(x, table, nomatch = 0L) == 0L


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### reading shapefiles ####
sf_use_s2(F)
all_amphibians = st_read("shapefiles/all_amphibians_dissolved.shp")
all_reptiles = st_read("shapefiles/GARD_17/Gard_1_7_ranges_fixex.shp", crs = 'EPSG:4326')%>%
  dplyr::rename(Binomial = binomial) 
all_mammals = readRDS("shapefiles/MAMMALS/MAMMALS_complete.shp")
all_birds0 = st_read("shapefiles/BIRDS/BIRDS_presence1_origin12_seasonal12.shp", crs = 'EPSG:4326') %>%  #fixed geometries %>% 
  dplyr::mutate(SCI_NAME = sci_name)

all_reptiles2 = all_reptiles %>% 
  dplyr::mutate(SCI_NAME = Binomial) 
all_mammals1 = all_mammals %>% 
  dplyr::rename("geometry" = "x")

all_birds_sp = all_birds0$sci_name %>% 
  as.data.frame()

df_list <- list(all_reptiles2 = all_reptiles2, all_amphibians = all_amphibians,all_mammals1= all_mammals1, all_birds0 = all_birds0)
all_vertebrates = bind_rows(df_list, .id = "source") %>% 
  dplyr::select(source, SCI_NAME) %>%
  mutate(source = recode(source, 
                         "all_amphibians" = "amphibians",
                         "all_reptiles2" = "reptiles",
                         "all_mammals1"= "mammals", 
                         "all_birds0" = "birds"))

continents = st_read("shapefiles/continents.shp")


### databases ####
LPI_populations = read_csv( "analysis/results/ALLvertebrates_dens_dist2_temporal_LPI_03032025.csv") %>% 
  dplyr::mutate(database = "LPI") 



library(dggridR)
library(collapse)
library(data.table)
dggs <- dgconstruct(res=12)

#reading the LPI dataframe after looping avg lambdas and distance to edges
LPI_populations$cel_id = dgGEO_to_SEQNUM(dggs, LPI_populations$lon, LPI_populations$lat)$seqnum

#reading the BIOTIME dataframe after looping avg lambdas and distance to edges
density_vertebrates_BIOTIME0 = read_csv("gridded_records_BIOTIME.csv")

#checking if there are duplicate populations (LPI points within BioTIME grid cells, with resolution defined by BioTIME)
common_species_grid = LPI_populations %>% 
  dplyr::right_join(as_data_frame(density_vertebrates_BIOTIME0 %>% 
                                    dplyr::mutate(Binomial = gsub("_", " ", Binomial))) , by = c("Binomial", "cel_id")) %>% 
  dplyr::filter(!is.na(database)) %>% 
  dplyr::mutate(geometry = as.character(geometry)) %>% 
  dplyr::select(Binomial, cel_id) %>% 
  dplyr::distinct()
# 2 species: Anaxyrus americanus and Larus americanus


BIOTME_populations = read_csv( "analysis/results/ALLvertebrates_dens_dist2_temporal_BIOTIME_10032025.csv") %>% 
  dplyr::left_join(density_vertebrates_BIOTIME0, by ="ID") %>% 
  dplyr::select(-"Binomial.y", -"Latitude", -"Longitude", -"count") %>% 
  dplyr::rename(Binomial= Binomial.x) %>% 
  dplyr::mutate(database = "BIOTIME") %>% 
  dplyr::anti_join(common_species_grid, by = c("Binomial", "cel_id")) %>%  # this eliminates duplicate records already identified
  dplyr::distinct(Binomial, cel_id, .keep_all = T)


hist(LPI_populations$avg_lambda, breaks = 100)
mean(LPI_populations$avg_lambda)

# creating a dataframe with both databases
database_populations = BIOTME_populations %>% 
  bind_rows(LPI_populations) 


database_populations%>% 
  distinct(Binomial, database) %>% 
  add_count(Binomial, name = "nr_databases") %>% 
  ggplot(., aes(x = nr_databases)) + geom_histogram() + theme_classic()

database_populations%>% 
  dplyr::distinct(Binomial, database) %>% 
  add_count(Binomial, name = "nr_databases") %>% 
  filter(nr_databases > 1) %>% 
  distinct(Binomial) %>% 
  dim()
#314 species common to both studies


### visualization + threats ####
library(raster)

database_populations_df = database_populations %>%
  mutate(increasing = ifelse(avg_lambda>0, 1, 0), # categorising population trends
         decreasing = ifelse(avg_lambda<0, 1, 0),
         lat_range = abs(ymin - ymax),
         dist_trailing = ifelse(lat>0, lat - ymin, ymax - lat), # calculating distance to trailing edges, corrected for latitude, CONFIRMED
         dist_leading = ifelse(lat>0, ymax - lat, lat - ymin), # calculating distance to leading edges, corrected for latitude, CONFIRMED
         reldist_leading = ifelse(lat>0, 
                                  ((lat - ymin) / (ymax - ymin)) * 100,  # calculating relative distance to leading edges, CONFIRMED 
                                  (abs(lat - ymax))/ (abs(ymin - ymax)) * 100), 
         reldist_trailing = 100 - reldist_leading, # calculates  relative distance to trailing edges as the inverse of reldist to leading 
         Group = case_when(Binomial %in% all_amphibians$SCI_NAME ~ "Amphibians",
                           Binomial %in% all_reptiles$Binomial~ "Reptiles",
                           Binomial %in% all_mammals$SCI_NAME~ "Mammals",
                           Binomial %in% all_birds0$sci_name~ "Birds"),
         Group = as.factor(Group)) %>% 
  add_count(Binomial, name = "nr_points") %>% 
  group_by(Binomial) %>%
  dplyr::mutate(correl = cor(avg_lambda,rgdc_poi, method = 'pearson', use = "na.or.complete"), # pearson correlation of distance to center and avg lambda, not used
                var_rgdcpoi = var(rgdc_poi)) %>% 
  dplyr::mutate(hemisphere_ymin = ifelse(ymin>0, "N", "S"),
                hemisphere_ymax = ifelse(ymax>0, "N", "S"),
                cross_equator = ifelse(hemisphere_ymin == hemisphere_ymax, 0, 1)) %>% # checking if species distribution crosses the equator
  dplyr::filter(lat_range > 5, #species with over 5º of latitudinal range, exclude range restricted species
                avg_lambda < 1.5,
                avg_lambda > -1.5) 

table(database_populations_df$Group)
database_populations_df %>% 
  # dplyr::filter(database == "BIOTIME") %>% 
  add_count(Binomial, name = "nr_points") %>%
  dplyr::filter(nr_points>2) %>%
  dplyr::distinct(Binomial, .keep_all = T) %>%
  # group_by(Group) %>% 
  # summarise(nr_spp = count(Group))
  dplyr::filter(Group == "Reptiles") %>% # excluding reptiles??
  View()

cols1 = c("Reptiles" = "olivedrab4",
          "Amphibians" = "turquoise4",
          "Mammals" = "chocolate4",
          "Birds" = "goldenrod1")
database_populations_df %>% 
  #  dplyr::filter(increasing == 0 & decreasing == 0) %>% 
  ggplot(., aes(x = reldist_trailing, fill = Group)) + geom_density(aes(alpha = 0.5))+ 
  scale_fill_manual(values = cols1) + theme_classic(base_size = 18)+
  labs(x = "Relative distance to trailing edge",
       y = "Density of records") + scale_alpha(guide = 'none')

database_populations_df %>% 
  dplyr::mutate(status = case_when(
    increasing == 0 & decreasing == 0 ~ "constant",
    increasing == 1 & decreasing == 0 ~ "increasing",
    increasing == 0 & decreasing == 1 ~ "decreasing")) %>% 
  # dplyr::filter(increasing == 0 & decreasing == 1) %>% 
  ggplot(., aes(x = reldist_trailing, fill = status)) + geom_density(alpha = 0.65)+ theme_classic(base_size = 18) +
  labs(x = "Relative distance to trailing edge",
       y = "Density of records") + facet_wrap(~Group) +
  scale_fill_manual(values = c("grey77", "darkred", "steelblue3"))

ggsave("preliminary_hist_reltrail.jpeg")

database_populations_df %>% 
  dplyr::filter(increasing == 1 & decreasing == 0) %>% 
  ggplot(., aes(x = reldist_leading, fill = Group)) + geom_density(aes(alpha = 0.5)) 

ggplot(database_populations_df, aes(reldist_trailing, avg_lambda)) + geom_point() +
  geom_smooth(method = "lm") + facet_wrap(~database) + theme_classic()

cols1 = c("Reptiles" = "olivedrab4",
          "Amphibians" = "turquoise4",
          "Mammals" = "chocolate4",
          "Birds" = "goldenrod1")

vertebrates_dens_temporal_populations2 = st_as_sf(database_populations_df, coords=c("lon", "lat"), crs= "EPSG:4326")
vertebrates_dens_temporal_populations3 =  st_intersects(vertebrates_dens_temporal_populations2, continents, sparse = F) #spatial join to get intersection of points and poly

populations_in_continents_vertebrates = vertebrates_dens_temporal_populations2[as.vector(vertebrates_dens_temporal_populations3), ]
populations_in_continents_vertebrates2 = vertebrates_dens_temporal_populations2 %>%  
  dplyr::mutate(Latitude = sf::st_coordinates(geometry)[,2],
                Longitude = sf::st_coordinates(geometry)[,1]) %>% 
  as.data.frame() %>% 
  dplyr::mutate(status = case_when(
    increasing == 0 & decreasing == 0 ~ "constant",
    increasing == 1 & decreasing == 0 ~ "increasing",
    increasing == 0 & decreasing == 1 ~ "decreasing"))

class(populations_in_continents_vertebrates2)

populations_in_continents_vertebrates %>% 
  dplyr::filter(avg_lambda<1.5,
                avg_lambda >-1.5) %>% 
  # dplyr::mutate(avg_lambda_10 = avg_lambda * 10) %>% 
  mapView(.)

write.csv(populations_in_continents_vertebrates2, "populations_in_continents_vertebrates.csv", row.names = FALSE)
populations_in_continents_vertebrates2 = read.csv("populations_in_continents_vertebrates.csv", row.names=NULL)

library("rnaturalearth")
world <- map_data("world")

#ideally change to a different crs
ggplot() +
  geom_polygon(data = world, aes(x=long, y=lat, group=group), colour="black",fill="grey", alpha=.3) + 
  theme_classic(base_size = 14) +
  geom_point(data=populations_in_continents_vertebrates2, aes(x=Longitude, y=Latitude, colour=status), pch=20, size=2.3, alpha = 0.6) +
  scale_color_manual(values = c("black", "darkred", "steelblue4")) +
  labs( x = "Longitude", y = "Latitude", colour = "Population trend") +
  coord_sf(ylim=c(-55,90)) + theme(legend.position = c(.15,.22))

ggsave("map_popstatus.jpeg")


example_decrease = populations_in_continents_vertebrates2 %>% 
  dplyr::mutate(lat_abs = abs(Latitude),
               lat_abs_scaled = scale(lat_abs),
               reldist_trailing_scaled = scale(reldist_trailing),
               reldist_leading_scaled = scale(reldist_leading)) %>% 
  dplyr::select(Binomial, Group, status, decreasing,
                reldist_trailing_scaled, lat_abs_scaled, last_year) %>% 
  slice_sample(n = 5000)

# write.csv(example_decrease, "example_decrease_Robin.csv")

head(example_decrease, n = 100)

example_decrease$status <- factor(example_decrease$status, levels = c("decreasing", "constant", "increasing"), ordered = TRUE)


ordered_model <- glmmTMB(
  formula = status ~ Group * (reldist_trailing_scaled + lat_abs_scaled) + (reldist_trailing_scaled | Binomial:Group) + (1 | last_year),
  data = example_decrease,
  family=ordbetareg # ordbetareg() is the family for ordered logit in glmmTMB
)
install.packages("glmmTMB",type="source")
library(glmmTMB)



decrease1 <- glmmTMB(formula = decreasing ~ reldist_trailing_scaled + 
                       (reldist_trailing_scaled|Binomial:Group) +
                       (1|last_year),
                     family = binomial,
                     data = example_decrease,
                     na.action = "na.fail")

summary(decrease1)
