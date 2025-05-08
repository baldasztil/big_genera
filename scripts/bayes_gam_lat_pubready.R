library(tidyverse)# fast csv reading
library(sf)
library(cowplot)
library(data.table)
library(ggpubr)

########### bayesian
library(brms)
library(tidybayes)
library(rstanarm)
library(marginaleffects)
library(broom)
library(broom.mixed) 
library(collapse)
library(WorldFlora)
library(patchwork)



# import datasets


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>%
  # remove countries without species
  filter(!LEVEL3_COD == "BOU") 


# native species 
dist_native <- fread("data/wcvp/dist_native.txt") 

# accepted names 
plants_full <- fread("data/wcvp/wcvp_accepted_merged.txt")


# big genera list 
big_genera <- read.csv("data/twenty_years_big.csv") 


midpoints_red <- read.table("data/midpoints_coordinates.txt") %>% 
  dplyr::select(LEVEL3_COD, lon, lat)

sp_info <- plants_full %>% 
  dplyr::select(plant_name_id, genus, family)

data(vascular.families)

angiosperms <- vascular.families %>% 
  filter(Group == "angiosperms")


# baseline richness ------------------------------------------------------------

# calculating overall richness patterns 
plants_big <- plants_full %>% 
  filter(genus %in% big_genera$Genus) 

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
  filter(family %in% angiosperms$Family)



# calculate richness patterns 
richness_patterns_bru <- dist_full %>% 
  filter(!Genus %in% c("Rubus", "Taraxacum", "Hieracium")) %>% 
  group_by(area_code_l3) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    big_sp = sum(presence > 0),
    non_big_sp = sum(presence < 0),
    prop_big = big_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) %>% 
  rename(LEVEL3_COD = area_code_l3)



percont_stats <- richness_patterns_bru %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  dplyr::select(-c(full_kcz, geometry)) %>% 
  as.data.frame() 




priors <- c(set_prior("normal(0.5, 0.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"), 
            set_prior("normal(0, 1)", class = "b", dpar = "phi"), 
            set_prior("normal(-4.6, 1)", class = "Intercept", dpar = "zi")
)


percont_stats$scaled_total_sp <- scale(percont_stats$total_sp)[,1]

bayes_model <- brm(bf(prop_big ~ s(lat), 
             phi ~ s(lat),
             zi ~ 1),
          data = percont_stats, 
          prior = priors, 
          cores = 4, seed = 17,
          family = zero_inflated_beta(),
          chains = 4, iter = 4000, warmup = 1000, 
          control = list(max_treedepth = 15, adapt_delta = 0.99))


options(scipen = 999)
full_slopes <- slopes(bayes_model)
hist(abs(full_slopes$estimate))
mean(abs(full_slopes$estimate))
sd(abs(full_slopes$estimate))
max(full_slopes$estimate)
min(full_slopes$estimate)

#mutate(across(where(is.numeric), ~ . * 100))

slopes_effect <- ggplot(full_slopes, aes(x = lat, y = estimate)) +
  geom_ribbon(full_slopes, mapping = aes(ymin = conf.low, ymax = conf.high), fill = "#cccccc",  alpha = 0.6) +
  geom_line(linewidth = 1.2, col = "red")+
  geom_hline(yintercept = 0, lty = 5, lwd = 0.8, col = "black") +
  labs(x = "Latitude", y = "Estimated effect") +
  theme_bw(base_size = 12) +
  coord_flip() 

  
bayes_extract_slopes <- bayes_model %>% 
    slopes(variables = "lat",
           newdata = datagrid(lat = c(-90, -60,-30, 0, 30, 60, 90))) %>% 
    posterior_draws()



effect_plot <- ggplot(bayes_extract_slopes, aes(x = draw, y = factor(lat), fill = factor(lat))) +
  geom_vline(xintercept = 0)  +
  stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_hdi",
               slab_alpha = 0.75, col = "black") +
  #scale_x_continuous(labels = label_pp_tiny) +
  scale_fill_viridis_d(option = "viridis", end = 1) +
  labs(x = "Average marginal effect on proportion of species in big plant genera", 
       y = "Latitude", fill = "Latitude",
       caption = "80% and 95% credible intervals shown in black") +
  theme_bw() +
  theme(legend.position = "bottom")



preds <- percont_stats %>%
  add_epred_draws(bayes_model)



bayes_gam <-
  ggplot(data = percont_stats, aes(x = lat, y = prop_big)) +
  stat_lineribbon(preds, mapping = aes(y = .epred), .width = c(.95, .8, .5),  col = "gold") +
  geom_point(aes(size = prop_big), alpha = 0.15, col = "#000c2f") +
  scale_fill_grey(start = 0.8, end = 0.4) +
  theme_bw(12) +
  coord_flip() +
  labs(x = "Latitude", y =  "Proportion of species in big plant genera", col = "Climate zones", fill = "Confidence level", 
       size = "Proportion")



model_plot <- bayes_gam + slopes_effect + plot_annotation(tag_levels = "A")

ggsave("cry_in_Lines_marginaleffect_GAM_bayes_prop_big_latitude.png", model_plot, width = 12, height = 5, dpi = 600)

model_plot
