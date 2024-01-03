library(tidyverse)# fast csv reading
library(sf)
library(rWCVP)
library(tmap)
library(tmaptools)


plotting.maps <- function (x) {
  
  mapping <- plants_megadiverse %>%
    filter(genus == x) 
  
  mapping_2 <- dist_megadiverse %>%  
    filter(plant_name_id %in% mapping$plant_name_id) %>% 
    group_by(area_code_l3) %>% 
    summarise(total_sp = n_distinct(plant_name_id)) %>%  
    rename(LEVEL3_COD = area_code_l3) 
  
  
  mapping_3 <- tdwg_3 %>% 
    left_join(mapping_2, by = "LEVEL3_COD") %>% 
    replace(is.na(.), 0)
  
  
  map <- tm_shape(mapping_3) + 
    tm_fill(col = "total_sp", 
            palette="Greys", 
            n = 5) +
    tm_legend(outside=TRUE) +
    tm_borders(col = "black", 
               lwd = 0.5, 
               lty = "solid",
               alpha = 0.3) +
    tm_layout(main.title  = paste0(x, " (n = ", length(unique(mapping$plant_name_id)), ")"),
              fontfamily = "serif", 
              main.title.size = 1.2, 
              main.title.position = c("left", "top")
    ) 
  
  tmap_save(map, paste0(x, "_map.svg"))
  
}

# this is the file with the species information
wcvp_raw <- read.table("data/wcvp/wcvp_names.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
# this file contains the distribution of each species  
dist_raw <- read.table("data/wcvp/wcvp_distribution.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 

mega_genera <- read.csv("data/mega_genera.csv")

# if you have species and distributions in one file it should also work! 
# Would just need to reformat some things / you could just generate a data frame
# that contains all the species names. 

# country codes 
tdwg_codes <- wgsrpd3 
tdwg_3 <- wgsrpd3 

# shapefile with midpoints of all countries 

# manipulate data --------------------------------------------------------------

# remove names without accepted name
wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

# remove all hybrids and genus names 
wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid =="")

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

# create combined dataset
plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")

# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 
plants_megadiverse <- plants_full %>% 
  filter(genus %in% mega_genera$megadiverse)

dist_megadiverse <- dist_native %>% 
  filter(plant_name_id %in% plants_megadiverse$plant_name_id) %>%  
  mutate(mega = "yes",
         presence = 1)

dist_normal <- dist_native %>% 
  filter(!plant_name_id %in% dist_megadiverse$plant_name_id) %>%  
  mutate(mega = "no", 
         presence = 0)

dist_full <- rbind(dist_megadiverse, dist_normal)

richness_patterns_bru <- dist_full %>% 
  group_by(area_code_l3) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    megadiverse_sp = sum(presence),
    prop_megadiverse = megadiverse_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) %>% 
  rename(LEVEL3_COD = area_code_l3)

richness_mapping <- tdwg_3 %>% 
  left_join(richness_patterns_bru, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0)



total_richness_map <- tm_shape(richness_mapping) + 
    tm_fill(col = "total_sp", 
            palette="Greens", 
            n = 5) +
    tm_legend(outside=TRUE) +
    tm_borders(col = "black", 
               lwd = 0.5, 
               lty = "solid",
               alpha = 0.3) +
    tm_layout(main.title  = "Total number of species",
              fontfamily = "serif", 
              main.title.size = 1.2, 
              main.title.position = c("left", "top")
    ) 

tmap_save(total_richness_map, "total_richness_map.svg")


megadiverse_richness_map <- tm_shape(richness_mapping) + 
  tm_fill(col = "megadiverse_sp", 
          palette="Oranges", 
          n = 5) +
  tm_legend(outside=TRUE) +
  tm_borders(col = "black", 
             lwd = 0.5, 
             lty = "solid",
             alpha = 0.3) +
  tm_layout(main.title  = "Total number of species in megadiverse genera",
           fontfamily = "serif", 
           main.title.size = 1.2, 
           main.title.position = c("left", "top")
           ) 

tmap_save(megadiverse_richness_map, "megadiverse_richness_map.svg")


prop_megadiverse_richness_map <- tm_shape(richness_mapping) + 
  tm_fill(col = "prop_megadiverse", 
          palette="Blues", 
          n = 5) +
  tm_legend(outside=TRUE) +
  tm_borders(col = "black", 
             lwd = 0.5, 
             lty = "solid",
             alpha = 0.3) +
  tm_layout(main.title  = "Proportion of species in megadiverse genera",
            fontfamily = "serif", 
            main.title.size = 1.2, 
            main.title.position = c("left", "top")
  ) 

tmap_save(prop_megadiverse_richness_map, "proportion_map.svg")


names <- mega_genera$megadiverse
names <- c("Begonia", "Solanum")

lapply(names, plotting.maps)
# Analysis ---------------------------------------------------------------------
summary(richness_patterns_bru)
boxplot(richness_mapping$prop_megadiverse)
boxplot(richness_mapping$megadiverse_sp)
plot(log10(richness_mapping$megadiverse_sp) ~ log10(richness_mapping$total_sp))
richness_mapping$log10megadiverse_sp <- log10(richness_mapping$megadiverse_sp)
richness_mapping$log10total_sp <- log10(richness_mapping$total_sp)

richness_stat_adj <- richness_mapping %>% 
  replace(is.na(.), 0) %>% 
  mutate(across(.cols = c("log10megadiverse_sp", "log10total_sp"), ~ ifelse(is.infinite(.x), 0, .x)))


a <-lm(richness_stat_adj$log10megadiverse_sp ~ richness_stat_adj$log10total_sp)
summary(a)
plot(a)
hist(log10(richness_stat_adj$megadiverse_sp))

