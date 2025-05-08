library(tidyverse)
library(sf)
library(rWCVP)
library(data.table)
library(WorldFlora)

# Baseline data ----------------------------------------------------------------

### To use this script a local copy of the WCVP data is needed. The data is openly available from: https://sftp.kew.org/pub/data-repositories/WCVP/ - wcvp v12 was used for the analysis here. 

# this is the file with the species information
wcvp_raw <- read.table("insert path here", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# this file contains the distribution of each species  
dist_raw <- read.table("insert path here", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 

# shapefile
tdwg_3 <-  st_read(dsn ="data/tdwg_3_realms_area/")  %>% 
  filter(!LEVEL3_COD == "BOU") %>% 
  rename(climate_zone = kcl_zone, realm = PhyloRealm) %>% 
  st_transform(crs = "+proj=eqearth")

# big genera list
big_genera <- read.csv("data/twenty_years_big_extract.csv")

# angiosperms
data("vascular.families")

angiosperms <- vascular.families %>% 
  filter(Group == "angiosperms")

# baseline richness ------------------------------------------------------------

# remove names without accepted name
wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

# remove all hybrids and genus names 
wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid =="") #%>% 
#filter(!genus %in% c("Hieracium", "Taraxacum", "Rubus", "Alchemilla"))

# filter distrbutions 
dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% # not introduced 
  filter(!location_doubtful == 1) %>%  # not doubtful occurrence 
  filter(!area == "") %>% # not unknown
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)

# summarise information
dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

# create combined angiosperm dataset
plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id") %>% 
  filter(family %in% angiosperms$Family) 

sp_info <- plants_full %>% 
  dplyr::select(plant_name_id, genus, family)

# calculating overall richness patterns 
plants_big <- plants_full %>% 
  filter(genus %in% big_genera$Genus) %>% 
  filter(family %in% angiosperms$Family)


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

fwrite(richness_mapping %>% st_drop_geometry(), "data/tdwg_overview_table_big_gen.csv")

size_big_genus <- dist_big %>% 
  group_by(genus) %>% 
  summarise(size =  n())


genus_matrix <- dist_big %>% 
  filter(!genus %in% c("Rubus", "Taraxacum", "Hieracium")) %>% 
  group_by(genus, realm) %>% 
  summarise(n = length(unique(plant_name_id))) %>% 
  left_join(size_big_genus, by = "genus") %>% 
  mutate(prop = n / size) %>% 
  dplyr::select(c(-n, -size)) %>% 
  pivot_wider(values_from = prop, 
              names_from = realm) %>% 
  replace(is.na(.),0) 

fwrite(genus_matrix, "data/big_genera_porportions_matrix.csv")


