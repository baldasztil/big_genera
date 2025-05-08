library(tidyverse)# fast csv reading
library(sf)
library(cowplot)
library(data.table)
library(rstatix)
library(ggpubr)
library(patchwork)
library(plotrix)
library(WorldFlora)
library(rnaturalearth)

# Data import ------------------------------------------------------------------

tdwg_3 <-  st_read(dsn ="data/tdwg_3_realms_area/")  %>% 
  st_transform(crs = "+proj=eqearth") %>% 
  dplyr::select(geometry, LEVEL3_COD)

richness_patterns_bru <- fread("data/tdwg_overview_table_big_gen.csv")

richness_mapping <- tdwg_3 %>% 
  left_join(richness_patterns_bru, by = "LEVEL3_COD")

# Analysis ---------------------------------------------------------------------


percont_stats <- richness_patterns_bru %>% 
  left_join(tdwg_3, by = c("LEVEL3_COD")) %>% 
  mutate(
         perc_big = prop_big * 100)  



dunn_prop_climate <- percont_stats %>% 
  dunn_test(prop_big ~ climate_zone) %>%   
  add_significance("p") %>% 
  add_y_position()



pal <- c( "#F1CE63", "#B07AA1", "#88CCEE",  "#EE8866","#2CA02C")


### Violin plots 
violins_prop <-  ggplot(percont_stats, aes(y = prop_big, x = reorder(climate_zone, lat))) +
  geom_jitter(width = 0.25, alpha = 0.5, aes(fill = climate_zone, col = climate_zone)) +
  geom_violin(alpha = 0.5, aes(fill = climate_zone, col = climate_zone)) +
  stat_pvalue_manual(dunn_prop_climate,label = "p.signif", tip.length = 0.01, size = 2.5) +
  theme_bw(12) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) + 
  labs(fill='Climate zone', col  = "Climate zone", x = "Climatic zone" , y = "Proportion of species in big genera") 
violins_prop

### Map of climate zones
sphere <- ne_download(category = "physical", type = "wgs84_bounding_box", returnclass = "sf")
sphere_trans <-  st_transform(sphere, st_crs(tdwg_3))

map_climate <-ggplot(richness_mapping) +
  geom_sf(data = sphere, fill = "white", col = "black", show.legend = F, lwd = 0.5) +
  geom_sf(aes(fill = climate_zone), col = "black", show.legend = F) +
  
  coord_sf(expand = T) +
  scale_fill_manual(values = pal) +
  labs(fill = "Climatic zone") + 
  theme_pubclean(base_size = 12)  



design <-   "
  11111
  11111
  11111
  22222
  22222
"


### Combined plot
(full_plot <- map_climate  +violins_prop + plot_layout(design = design, heights = c(0.3,0.3,0.3,1,0.1), guides = "collect")  + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold')))

