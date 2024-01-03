library(tidyverse)# fast csv reading
library(sf)
library(rWCVP)
library(tmap)
library(tmaptools)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
library(viridis)
library(FSA)
library(plotrix)
library(ggcorrplot)
library(treemap)
library(treemapify)




#abce

# this is the file with the species information
wcvp_raw <- read.table("data/wcvp/wcvp_names.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
# this file contains the distribution of each species  
dist_raw <- read.table("data/wcvp/wcvp_distribution.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 

mega_genera <- read.csv("data/mega_genera.csv")

midpoints_red <- read.table("data/midpoints_coordinates.txt") %>% 
  dplyr::select(LEVEL3_COD, lon, lat)

tdwg_codes <- read.table("data/tdwg_codes.csv", sep = ",", header = T) %>% 
  dplyr::select(LEVEL3_COD, LEVEL1_NAM, LEVEL1_COD) %>% 
  st_drop_geometry()


# if you have species and distributions in one file it should also work! 
# Would just need to reformat some things / you could just generate a data frame
# that contains all the species names. 

# country codes 
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
  #filter(!genus %in% c("Hieracium", "Taraxacum"))

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

# score megadiverse genera 
dist_megadiverse <- dist_native %>% 
  filter(plant_name_id %in% plants_megadiverse$plant_name_id) %>%  
  mutate(mega = "yes",
         presence = 1)

# score non-megadiverse genera
dist_normal <- dist_native %>% 
  filter(!plant_name_id %in% dist_megadiverse$plant_name_id) %>%  
  mutate(mega = "no", 
         presence = -1)

# combine dataset
dist_full <- rbind(dist_megadiverse, dist_normal) %>%
  left_join(wcvp_accepted, by = "plant_name_id")
  

# calculate richness patterns 
richness_patterns_bru <- dist_full %>% 
  group_by(area_code_l3) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    megadiverse_sp = sum(presence > 0),
    non_megadiverse_sp = sum(presence < 0),
    prop_megadiverse = megadiverse_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) %>% 
  rename(LEVEL3_COD = area_code_l3)


# reformat for mapping
richness_mapping <- tdwg_3 %>% 
  left_join(richness_patterns_bru, by = "LEVEL3_COD") %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") %>% 
  left_join(tdwg_codes, by = "LEVEL3_COD") %>% 
  replace(is.na(.), 0) 


# plotting richness patterns ----------------------------------------------------

# visualising the latitudinal gradient with scatterplots / bubbleplots

a <- ggplot(richness_mapping, aes(x = total_sp, y = lat)) +
  geom_point(color = "darkgreen", alpha = 0.3) +
  labs(x = "Number of species", y = "Latitude") +
  geom_smooth(orientation = "y", col = "black", se = T, alpha = 0.15) +
  theme_bw()

b <- ggplot(richness_mapping, aes(x = megadiverse_sp, y = lat)) +
  geom_point(color = "darkorange", alpha = 0.3) +
  labs(x = "Number of species in megadiverse genera", y = "Latitude") +
  geom_smooth(orientation = "y", col = "black", se = T, alpha = 0.15)+
  theme_bw() 

c <- ggplot(richness_mapping, aes(x = prop_megadiverse, y = lat)) +
  geom_point(color = "darkblue", alpha = 0.3) +
  geom_smooth(orientation = "y", col = "black", se = T, alpha = 0.15) +
  labs(x = "Proportion of species in megadiverse genera", y = "Latitude") +
  theme_bw() 

# correlation between megadiversity and diversity with proportion as size
d <- ggplot(richness_mapping, aes(x = megadiverse_sp, y = total_sp)) +
  geom_point(aes(size = prop_megadiverse, col = LEVEL1_NAM),  alpha = 0.6) +
  geom_text_repel(aes(label = LEVEL3_COD),
                   segment.color ='grey',  size = 2.5) +
  scale_color_brewer("Continent", palette = "Set3") +
  geom_smooth(orientation = "x", col = "black", se = T, alpha = 0.15) +
  labs(x = "Number of species in megadiverse genera", y = "Number of species") +
  scale_size("Proportion in mega genera") +
  theme_bw() 

# layout 
e <- a + b + c + plot_layout(ncol = 1, nrow = 3)

# layout 2
f <- e | d +  plot_layout(guides = "collect")

ggsave(file = paste0("latidunal_patterns.svg"),f,   
       width = 12, height = 7)



# Testing the correlation between mega diversity and continent -----------------

# are patterns  significantly different between continents?
a <- dunnTest(prop_megadiverse ~ LEVEL1_NAM,
         data = richness_mapping,
         method = "bonferroni", 
         two.sided = T)$res

# reformating data 
wide_format <-  a %>% 
  dplyr::select(Comparison, P.adj) %>% 
  mutate(Comparisons = Comparison) %>% 
  separate(Comparison, into = c("Continent1", "Continent2"), sep = " - ")

# continued with loop to transform results into correlation matrix 
c <- unique(richness_mapping$LEVEL1_NAM)
c <- sort(c)

resultlist <- list() 

for (i in 1:9) {
  continent <- c[i]
  
  temp <- wide_format %>% 
    filter(Continent1 == continent | Continent2 == continent) %>% 
    mutate(Continent3 = case_when(!Continent1 == continent~ Continent1, 
                                TRUE ~ Continent2),  
           Continent1 = continent) %>% 
   #mutate(P.adj = case_when(P.adj < 0.05 ~ 1, 
    #                         TRUE ~ 0)) %>% 
    dplyr::select(Continent1, Continent3, P.adj) %>% 
    
    pivot_wider(values_from = P.adj, names_from = Continent3) 
  
  resultlist[[i]] <- temp
}

results_lon <- as.data.frame(do.call(bind_rows, resultlist))
rownames(results_lon) <- results_lon$Continent

cor_matrix <- results_lon %>%
  select(order(colnames(.))) %>% 
  dplyr::select(-Continent1) %>% 
  as.matrix() 

# visualisation as correlation plot
corrplot <-ggcorrplot(cor_matrix, p.mat = cor_matrix, insig = "pch", method = "square", tl.cex = 9, 
           title = "Significance of correlation between continent 
           and proportion of megadiverse sp",
           legend.title = "p-value", show.diag = T, hc.order = F) + 
  scale_fill_viridis(limits = c(0, 0.05), na.value = "white") 

# creating a heatmapt of the patterns  -----------------------------------------
# 
richness_mapping_long <- richness_mapping %>% 
  dplyr::select(total_sp, megadiverse_sp, prop_megadiverse, LEVEL1_NAM, LEVEL3_COD, lat) %>% 
  mutate(across(c(total_sp, megadiverse_sp, prop_megadiverse), scale)) %>% 
  pivot_longer(cols = c(total_sp, megadiverse_sp, prop_megadiverse), 
               names_to = "Index", 
               values_to = "Value") %>% 
  st_drop_geometry() 

# reordering the data 
richness_mapping_long$LEVEL3_COD <- reorder(richness_mapping_long$LEVEL3_COD, 
                                            richness_mapping_long$lat)

# visualisation 
heatmap <- ggplot(richness_mapping_long, aes(x = LEVEL3_COD, y = Index, fill = Value, col = Value)) +
  geom_tile(height = 1.5, width =1.1,  linejoin = "round", linewidth = 0, color = NA) + 
  scale_fill_viridis_c() + 
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 4)) +
  xlab("Country code") +
  ylab("Scaled mean") + 
  ggtitle("Latitudinal patterns of diversity") +
  coord_flip()


ggsave("heatmap.svg", width = 20, 
       height = 20)

# continental patterns  --------------------------------------------------------

# summarising the richness patterns 
continent_patterns <- as.data.frame(richness_mapping) %>%
  group_by(LEVEL1_NAM) %>% 
  summarise(
            mean_sp = mean(total_sp), 
            sd_total = std.error(total_sp),
            
            mean_mega = mean(megadiverse_sp), 
            sd_mega = std.error(megadiverse_sp),
            
            mean_nonmega = mean(non_megadiverse_sp), 
            sd_nonmega = std.error(non_megadiverse_sp),
      
            mean_prop = mean(prop_megadiverse), 
            sd_propmega = std.error(prop_megadiverse),
            
            count = n()
            )

# reformatting the data 
continent_patterns_long <- continent_patterns %>% 
  pivot_longer(cols = c("mean_mega","mean_nonmega", "sd_mega", "sd_nonmega"), 
                 names_to = c(".value", "megadiverse"),
               names_sep = "_") %>% 
  dplyr::select(Continent = LEVEL1_NAM, group = megadiverse, Species = mean, sd)


# visualisation as barplot 
barplot <- ggplot(continent_patterns_long, aes(fill=group, y= Continent, x= Species)) + 
  geom_bar(position="dodge", stat="identity", alpha = 0.9) +
  scale_fill_viridis(discrete = T, begin = 1, end = 0) +
  geom_errorbar(aes(y = Continent, xmin= Species-sd, xmax= Species+sd), width=0.5, 
                colour="black", alpha=0.95, size=0.75, 
                position = position_dodge(0.9)) +
  theme_bw() + 
  xlab("Mean number of species") +
  ggtitle("Species means per continent")

ggsave("barplot.svg", height = 10, width = 10)


# combining the elements 
bar_cor <- barplot / corrplot + plot_layout(nrow = 2)
heatbarcor <- bar_cor | heatmap 

ggsave("heatbarcor.svg",heatbarcor, width = 20, 
       height = 20)


# treemap trials  --------------------------------------------------------------

treemapping <- richness_mapping %>% 
  dplyr::select(LEVEL3_COD, LEVEL1_NAM, total_sp, megadiverse_sp, non_megadiverse_sp, 
                prop_megadiverse) %>% 
  st_drop_geometry()
  
treemap(treemapping,
        index=c("LEVEL1_NAM", "LEVEL3_COD"),
        vSize="total_sp",
        type="index"
)

treemap(treemapping,
        index=c("LEVEL1_NAM", "LEVEL3_COD"),
        vSize="megadiverse_sp",
        type="index"
)
svg("treemap_nogg.svg", width = 10, height = 10)
treemap(treemapping,
        index=c("LEVEL1_NAM", "LEVEL3_COD"),
        vSize="prop_megadiverse",
        type="index", 
        algorithm = "pivotSize",
        fontsize.labels=c(12,8),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","grey15"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("center", "center")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=1,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
        
        title  = "Proportion of species from megadiverse genera", 
        palette = "Set3", 
        border.lwds	= NA, 
        border.col = "white")
dev.off()

library(ggplot2)
library(treemap)
tree_map <- ggplot(treemapping, aes(area = prop_megadiverse, fill = LEVEL1_NAM, label = LEVEL3_COD,
                     subgroup = LEVEL1_NAM)) +
  geom_treemap(show.legend = F, col = "white")+
  geom_treemap_text(grow = F, reflow = F, size = 9, place = "centre", alpha = 1) + 
  geom_treemap_subgroup_border(show.legend = F, col = "black", size = 3) +
  geom_treemap_subgroup_text(place = "centre", grow = F,
                             reflow = F, 
                             alpha = 0.7, colour = "black",
                             fontface = "bold", 
                             size = 14) +
  scale_fill_brewer(palette = "Set3") + 
  ggtitle("Proportion of species in megadiverse genera")

ggsave("treemap.svg",tree_map, width = 10, height = 10)






# circular bar plot of top 5 genera per continent-------------------------------
continent_genus_patterns <- dist_full %>% 
  group_by(continent, genus, mega) %>% 
  summarise (sp = n_distinct(plant_name_id)) %>% 
  group_by(continent) %>% 
  slice_max(n = 5, order_by = sp)




# Create dataset
data <- continent_genus_patterns %>% 
  dplyr::select(individual = genus, group = continent, value = sp) %>% 
  mutate(group = as.factor(group))

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
circular_barplot <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 500, xend = start, yend = 500), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1000, xend = start, yend = 1000), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1500, xend = start, yend = 1500), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(100, 500, 1000, 1500), label = c("100", "500", "1000", "1500") , color="grey", size=3 , angle=0, fontface="bold", hjust=1.5) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-2000,max(data$value)+10) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface=c("bold.italic"),alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start - 0.5, y = -30, xend = end + 0.5, yend = -30), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -470, label=group), hjust=c(0.5,0.5,0.5,0.5,0.5,1,0.4,0.4,0.5), colour = "black", alpha=0.8, size= 2, fontface="bold", inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("The five biggest genera per continent")


ggsave(file = paste0("circulat_barplot.svg"),circular_barplot,   
       width = 20, height = 18)






# combining plots --------------------------------------------------------------



prop_analysis <- tree_map + circular_barplot + plot_layout(guides = "keep")

ggsave(file = paste0("latitudinal_patterns_2.svg"),prop_analysis,   
       width = 20, height = 20)





# test -------------------------------------------------------------------------



richness_patterns_test <- dist_full %>% 
  filter(area_code_l3 == "SWE") %>% 
  group_by(mega) %>% 
  summarise(
    total_sp = n_distinct(plant_name_id),
    megadiverse_sp = sum(presence),
    prop_megadiverse = megadiverse_sp / total_sp
  ) %>%
  arrange(desc(total_sp)) 


dist_test <- dist_full %>% 
  filter(area_code_l3 == "SWE") 


plants_test <- plants_full %>% 
  filter(plant_name_id %in% dist_test$plant_name_id) %>% 
  group_by(genus) %>% 
  summarise(n = length(unique(plant_name_id)))
