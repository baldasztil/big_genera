library(tidyverse)# fast csv reading
library(sf)
library(cowplot)
library(data.table)
library(ggpubr)
library(brms)
library(tidybayes)
library(rstanarm)
library(marginaleffects)
library(broom)
library(broom.mixed) 
library(collapse)
library(WorldFlora)
library(patchwork)



# Data import ------------------------------------------------------------------


tdwg_3 <-  st_read(dsn ="data/tdwg_3_realms_area/")  %>% 
  filter(!LEVEL3_COD == "BOU") %>% 
  rename(climate_zone = kcl_zone, realm = PhyloRealm) %>% 
  st_transform(crs = "+proj=eqearth")

richness_patterns_bru <- fread("data/tdwg_overview_table_big_gen.csv")


# Modelling --------------------------------------------------------------------


priors <- c(set_prior("normal(0.5, 0.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"), 
            set_prior("normal(0, 1)", class = "b", dpar = "phi"), 
            set_prior("normal(-4.6, 1)", class = "Intercept", dpar = "zi")
)


richness_patterns_bru$scaled_total_sp <- scale(richness_patterns_bru$total_sp)[,1]

bayes_model <- brm(bf(prop_big ~ s(lat), 
             phi ~ s(lat),
             zi ~ 1),
          data = richness_patterns_bru, 
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


# Plots ------------------------------------------------------------------------ 

# Slopes and marginal effect 
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



# Smooth plot 

preds <- richness_patterns_bru %>%
  add_epred_draws(bayes_model)



bayes_gam <-
  ggplot(data = richness_patterns_bru, aes(x = lat, y = prop_big)) +
  stat_lineribbon(preds, mapping = aes(y = .epred), .width = c(.95, .8, .5),  col = "gold") +
  geom_point(aes(size = prop_big), alpha = 0.15, col = "#000c2f") +
  scale_fill_grey(start = 0.8, end = 0.4) +
  theme_bw(12) +
  coord_flip() +
  labs(x = "Latitude", y =  "Proportion of species in big plant genera", col = "Climate zones", fill = "Confidence level", 
       size = "Proportion")



(model_plot <- bayes_gam + slopes_effect + plot_annotation(tag_levels = "A")) 



