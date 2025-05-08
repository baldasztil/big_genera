library(tidyverse)# fast csv reading
library(sf)
library(tmap)
library(tmaptools)
library(data.table)
library(WorldFlora)



normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

# data import ------------------------------------------------------------------


tdwg_3 <-  st_read(dsn ="data/tdwg_3_realms_area/")  %>% 
  st_transform(crs = "+proj=eqearth") %>% 
  dplyr::select(geometry, LEVEL3_COD)

richness_patterns_bru <- fread("data/tdwg_overview_table_big_gen.csv")

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










  

