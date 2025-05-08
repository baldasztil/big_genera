library(tidyverse)# fast csv reading
library(sf)
library(tmap)
library(tmaptools)
library(data.table)
library(WorldFlora)



normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}


tdwg_3 <-  st_read(dsn ="data/wgsrpd-master/tdwg_3_realms_area/")  %>% 
  filter(!LEVEL3_COD == "BOU") %>% 
  rename(climate_zone = kcl_zone, realm = PhyloRealm) %>% 
  st_transform(crs = "+proj=eqearth")



dist_native <- fread("data/wcvp/dist_native.txt")  
big_genera <- read.csv("data/twenty_years_big.csv")
data("vascular.families")
  
angiosperms <- vascular.families %>% 
  filter(Group == "angiosperms")

plants_full <- fread("data/wcvp/wcvp_accepted_merged.txt") %>% 
  filter(family %in% angiosperms$Family) 

sp_info <- plants_full %>% 
  dplyr::select(plant_name_id, genus, family)
  

# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 
# calculating overall richness patterns 
plants_big <- plants_full %>% 
  filter(genus %in% big_genera$Genus) %>% 
  filter(family %in% angiosperms$Family)# %>% 
  #filter(!genus %in% c("Rubus", "Taraxacum", "Hieracium")) 


# score big genera 
dist_big <- dist_native %>% 
  filter(plant_name_id %in% plants_big$plant_name_id) %>%  
  mutate(big = "yes",
         presence = 1)

# score non-big genera
dist_normal <- dist_native %>% 
  filter(!plant_name_id %in% dist_big$plant_name_id) %>%  
  mutate(big = "no", 
         presence = -1)

# combine dataset
dist_full <- rbind(dist_big, dist_normal) %>%
  left_join(sp_info, by = "plant_name_id")  %>% 
  filter(family %in% angiosperms$Family) %>% 
 filter(!genus %in% c("Rubus", "Taraxacum", "Hieracium")) 
  
length(unique(dist_full$plant_name_id))



# calculate richness patterns 
richness_patterns_bru <- dist_full %>% 
  group_by(area_code_l3) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    big_sp = sum(presence > 0),
    non_big_sp = sum(presence < 0),
    prop_big = big_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) %>% 
  rename(LEVEL3_COD = area_code_l3)

richness_mapping <- tdwg_3 %>% 
  left_join(richness_patterns_bru, by = "LEVEL3_COD")



# Analysis ---------------------------------------------------------------------



richness_mapping_long <- richness_mapping %>% 
  pivot_longer(names_to = "dataset", 
               values_to = "diversity", 
               cols = c("total_sp", "big_sp", "non_big_sp", "prop_big"))  %>% 
  mutate(dataset = factor(dataset, levels = c("total_sp", "non_big_sp", "big_sp", "prop_big"))) %>% 
  group_by(dataset) %>% 
  mutate(div_norm = normalize(diversity))

  
col_names <- c("Total species", 
               "Species in genera < 500 spp.", 
               "Species big plant genera", 
               "Proportion of species big plant genera")

(multiple_diversity_maps <- tm_shape(richness_mapping_long) + 
  tm_fill(col = "diversity", 
          palette="YlOrRd", 
          n = 5, style = "pretty") +
  #tm_format("World") +
  tm_borders() +
  tm_style("white", earth.boundary = c(-180, -87, 180, 87), 
           space.color = "white", earth.boundary.lwd = 2.5, legend.outside = F) +
  tm_facets(by = "dataset", free.scales = T, ncol = 2) +
  tm_layout(frame = F, frame.lwd = 0, panel.label.bg.color = "black", asp = 1.2, 
            panel.label.size = 1.5,  panel.label.height = 1,
            panel.labels = col_names, 
            panel.label.color = "white",
            legend.outside = F, 
            
  ))










  

