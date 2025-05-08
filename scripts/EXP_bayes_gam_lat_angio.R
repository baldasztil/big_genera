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


plot_post_linPred <- function(fit, data) {
  postPred <- data |> 
    tidybayes::add_epred_draws(fit) |> 
    group_by(prop_big,lat) |> 
    summarize(
      lower  = tidybayes::hdi(.epred)[1],
      mean   = mean(.epred),
      higher = tidybayes::hdi(.epred)[2]
    )
  postPred
}


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU") 

tdwg_3_area <- tdwg_3 %>% 
  mutate(area = as.numeric(st_area(.) / 1000000)) %>% 
  st_drop_geometry()

tdwg_codes <- tdwg_3 %>%  
  dplyr::select(LEVEL3_COD, LEVEL1_NAM, LEVEL1_COD, climate_zone = kcl_zone, full_kcz) %>% 
  st_drop_geometry()

dist_native <- fread("data/wcvp/dist_native.txt") 
plants_full <- fread("data/wcvp/wcvp_accepted_merged.txt")



big_genera <- read.csv("data/twenty_years_big.csv") 
# %>% 
#   filter(!Genus %in% c("Rubus", "Taraxacum", "Hieracium"))

midpoints_red <- read.table("data/midpoints_coordinates.txt") %>% 
  dplyr::select(LEVEL3_COD, lon, lat)

sp_info <- plants_full %>% 
  dplyr::select(plant_name_id, genus, family)

library(WorldFlora)
data(vascular.families)

angiosperms <- vascular.families %>% 
  filter(Group == "angiosperms")


# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 
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

length(unique(dist_full$plant_name_id))
# %>%  
#   filter(!genus %in% c("Rubus", "Taraxacum", "Hieracium"))


# calculate richness patterns 
richness_patterns_bru <- dist_full %>% 
  #filter(!genus %in% c("Rubus", "Taraxacum", "Hieracium")) %>% 
  group_by(area_code_l3) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    big_sp = sum(presence > 0),
    non_big_sp = sum(presence < 0),
    prop_big = big_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) %>% 
  rename(LEVEL3_COD = area_code_l3) 
  #mutate(prop_big = ifelse(prop_big == 0, 0.0000001, prop_big))
#filter(prop_big > 0)


percont_stats <- richness_patterns_bru %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  dplyr::select(-c(full_kcz, geometry)) %>% 
  as.data.frame() 

ggplot(percont_stats, aes( x = lat, y = prop_big)) +
  geom_point() +
  geom_smooth()



get_prior(
  prop_big ~ s(lat),
  data = percont_stats,
  family = zero_inflated_beta()
)



priors <- c(
  # Priors for the mean of the Beta distribution (proportion data)
  set_prior("normal(0, 1)", class = "b"),        # Fixed effects for predictors (e.g., latitude)
  set_prior("normal(0, 1)", class = "b", dpar = "phi"), 
  set_prior("normal(-1.3, 2)", class = "Intercept"),# Intercept prior
  # Priors for zero-inflation component (probability of zeroes)
  set_prior("normal(-4, 1)", class = "Intercept", dpar = "zi")  # Zero-inflation logit (probability of zeroes)
)


percont_stats$scaled_total_sp <- scale(percont_stats$total_sp)[,1]



bayes_model <- brm(bf(prop_big ~ s(lat), 
             phi ~ s(lat),
             zi ~ 1),
          data = percont_stats, 
          prior = priors, 
          cores = 4, seed = 20,
          family = zero_inflated_beta(),
          chains = 4, iter = 4000, warmup = 1000, 
          control = list( adapt_delta = 0.99))


default_priors <- c(
  prior(normal(0, 5), class = "b") # Prior for the residual standard deviation (sigma)
)


priors <- c(set_prior("normal(0.5, 0.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"), 
            set_prior("normal(0, 1)",  class = "sigma"), 
            set_prior("beta(0.1, 10)",   class = "Intercept", dpar = "hu")
)

priors <- c(set_prior("normal(0.5, 0.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"), 
            set_prior("normal(0, 1)",  class = "sigma"), 
            set_prior("beta(0.1, 10)",   class = "Intercept", dpar = "hu")
)


priors <- c(
  # Hurdle (zero-inflation) component: prior for logit of the probability of >0
  prior(normal(0, 1), class = "b", dpar = "hu"),
  
  # Log-normal part: priors for the log-transformed mean and variance
  prior(normal(0, 2),  dpar = "mu"),  # Log-scale mean
  prior(normal(0, 1),  dpar = "sigma"))

bayes_model <- brm(bf((prop_big) ~ s(lat)),
                   hu ~ 1,
                   data = percont_stats, 
                   cores = 4, seed = 17,
                   prior = default_priors, 
                   family = hurdle_gamma(),
                   chains = 4, iter = 4000, warmup = 1000, 
                   control = list(max_treedepth = 15, adapt_delta = 0.99))


summary(bayes_model)
plot(bayes_model)

bayes_model$prior

head(as_draws_df(bayes_model))



bayes_R2(bayes_model)

x <- as_draws_df(bayes_model)


percont_stats_red <- percont_stats %>% 
  filter(prop_big > 0)
fitdistrplus::descdist(percont_stats$prop_big)
x <- fitdistrplus::fitdist(log10(percont_stats$prop_big), "lnorm")
x <- fitdistrplus::fitdist(percont_stats$prop_big , "gamma")


plot(x)
hist(percont_stats$prop_big)

plot(density(rbeta(1000, 0.1,10)))

plot(density(rnorm(1000, mean= -1.3,sd= 2)))


pp_check(bayes_model, ndraws = 1000)
pp_check(bayes_model, type = "ecdf_overlay")
pp_check(bayes_model,
         type = "stat")


smooth_plot <- conditional_smooths(bayes_model)
smooth_plot

x <- posterior_epred(bayes_model)




beta_bayes_pred <- bayes_model %>% 
  epred_draws(newdata = tibble(lat = c(-84.836909,  -3.724873,  20.248999,  41.085648,  78.770896)), ndraws = 1000) %>% 
  mutate(lat = as.factor(lat))


ggplot(beta_bayes_pred, aes(x = .epred, y = lat, fill = lat)) +
  stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_hdi") +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  guides(fill = "none") +
  theme_bw()

beta_bayes_pred_1 <- bayes_model %>% 
  epred_draws(newdata = expand_grid(
                                    lat = seq(-90, 90, by = 1)))


ggplot(beta_bayes_pred_1, aes(x = lat, y = .epred)) +
  stat_lineribbon() + 
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Latitude", y = "Predicted proportion of species in genera > 500 spp.",
       fill = "Credible interval") +
  theme_bw() +
  theme(legend.position = "bottom")




full_slopes <- slopes(bayes_model)
hist(abs(full_slopes$estimate))
mean(abs(full_slopes$estimate))
sd(abs(full_slopes$estimate))
#mutate(across(where(is.numeric), ~ . * 100))

slopes_effect <- ggplot(full_slopes, aes(x = lat, y = estimate)) +
  geom_ribbon(full_slopes, mapping = aes(ymin = conf.low, ymax = conf.high), fill = "#cccccc",  alpha = 0.6) +
  geom_line(linewidth = 1.2, col = "red")+
  geom_hline(yintercept = 0, lty = 5, lwd = 0.8, col = "black") +
  labs(x = "Latitude", y = "Estimated effect") +
  theme_bw(base_size = 12) +
  coord_flip() 

slopes_effect
ggsave("cry_rem_marginal_effect_bayes_prop_big_latitude_angio.svg", slopes_effect, width = 7, height = 5)


  
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
  labs(x = "Average marginal effect on proportion of species in genera > 500 spp.", 
       y = "Latitude", fill = "Latitude",
       caption = "80% and 95% credible intervals shown in black") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("cry_rem_Density_marginal_effect_bayes_prop_big_latitude_angio.svg", effect_plot, width = 7, height = 5)




tidybayes::summarise_draws(bayes_model)



preds1 <- plot_post_linPred(bayes_model, percont_stats)


# ggplot(preds1, aes(x = prop_big, y = mean)) +
#   geom_point() +
#   coord_fixed() +
#   scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,0.8)) +
#   scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,0.8))
 
#cor.test(preds1$prop_big, preds1$mean)

preds <- percont_stats %>%
  add_epred_draws(bayes_model)

plot(percont_stats$prop_big ~ percont_stats$area)


pal <- nationalparkcolors::park_palette("Acadia")[c(1,2,3,5,6)]

bayes_gam <-
  ggplot(data = percont_stats, aes(x = lat, y = prop_big)) +
  stat_lineribbon(preds, mapping = aes(y = .epred), .width = c(.95, .8, .5),  col = "gold") +
  geom_point(aes(size = total_sp), alpha = 0.15, col = "#000c2f") +
  scale_fill_grey(start = 0.8, end = 0.4) +
  scale_color_manual(values = pal) +
  theme_bw(12) +
  coord_flip() +
  labs(x = "Latitude", y =  "Proportion of species in genera > 500 spp.", col = "Climate zones", fill = "Confidence level", 
       size = "Total species")

bayes_gam
ggsave("cry_rem_GAM_bayes_prop_big_latitude_angio.svg", bayes_gam, width = 7, height = 5)


library(patchwork)
model_plot <- bayes_gam + slopes_effect + plot_annotation(tag_levels = "A")

ggsave("cry_rem_Lines_marginaleffect_GAM_bayes_prop_big_latitude_angio.png", model_plot, width = 12, height = 5, dpi = 600)

model_plot
