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

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/tdwg_3_realms/") %>% 
  filter(!LEVEL3_COD == "BOU") %>% 
  rename(climate_zone = kcl_zone) %>%
  st_transform(crs = "+proj=eqearth") %>% 
  mutate(climate_zone = ifelse(climate_zone == "Cold", "Continental", climate_zone))

tdwg_3_area <- tdwg_3 %>% 
  mutate(area = as.numeric(st_area(.) / 1000000)) %>% 
  st_drop_geometry()

tdwg_codes <- tdwg_3 %>%  
  dplyr::select(LEVEL3_COD, LEVEL1_NAM, LEVEL1_COD, climate_zone, full_kcz) %>% 
  st_drop_geometry()

dist_native <- fread("data/wcvp/dist_native.txt") 
plants_full <- fread("data/wcvp/wcvp_accepted_merged.txt")



big_genera <- read.csv("data/twenty_years_big_extract.csv") #%>% 
 # filter(!Genus %in% c("Rubus", "Taraxacum", "Hieracium"))


midpoints_red <- read.table("data/midpoints_coordinates.txt") %>% 
  dplyr::select(LEVEL3_COD, lon, lat)

sp_info <- plants_full %>% 
  dplyr::select(plant_name_id, genus, family)

data(vascular.families)

angiosperms <- vascular.families %>% 
  filter(Group == "angiosperms")

# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 
# calculating overall richness patterns 
plants_big <- plants_full %>% 
  filter(genus %in% big_genera$Genus) 

length(unique(plants_big$genus))

# score big genera 
dist_big <- dist_native %>% 
  filter(plant_name_id %in% plants_big$plant_name_id) %>%  
  mutate(mega = "yes",
         presence = 1)

xx <- plants_big %>% 
  group_by(genus) %>% 
  summarise(n = n_distinct(plant_name_id))

# score non-big genera
dist_normal <- dist_native %>% 
  filter(!plant_name_id %in% dist_big$plant_name_id) %>%  
  mutate(mega = "no", 
         presence = -1)

# combine dataset
dist_full <- rbind(dist_big, dist_normal) %>%
  left_join(sp_info, by = "plant_name_id") %>% 
  filter(family %in% angiosperms$Family)


# calculate richness patterns 
richness_patterns_bru <- dist_full %>% 
  filter(!genus %in% c("Rubus", "Taraxacum", "Hieracium")) %>% 
  group_by(area_code_l3) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    big_sp = sum(presence > 0),
    non_big_sp = sum(presence < 0),
    prop_big = big_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) %>% 
  rename(LEVEL3_COD = area_code_l3)
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

map_climate <-ggplot(tdwg_3) +
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

