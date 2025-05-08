library(tidyverse)# fast csv reading
library(sf)
library(tmap)
library(tmaptools)
library(patchwork)
library(cluster)
library(factoextra)
library(data.table)
library(pvclust)
library(NbClust)
library(pheatmap)
library(viridis)
library(ggtree)
library(ape)


replace.cluster <- function(x) {
  
  x_cluster_na <- unique(x$cluster)
  x_cluster <- x_cluster_na[which(!is.na(x_cluster_na))]
  
  replacement <- x %>% 
    mutate(cluster = replace_na(cluster, x_cluster)) 
  
  return(replacement)
  
}


tdwg_3 <-  st_read(dsn ="data/tdwg_3_realms_area/")  %>% 
  st_transform(crs = "+proj=eqearth") %>% 
  dplyr::select(geometry, LEVEL3_COD)

richness_patterns_bru <- fread("data/tdwg_overview_table_big_gen.csv")

richness_mapping <- tdwg_3 %>% 
  left_join(richness_patterns_bru, by = "LEVEL3_COD")


genus_matrix <- fread("data/big_genera_porportions_matrix_pubready.csv")



# dataformatting ---------------------------------------------------------------

df <- genus_matrix %>% 
  ungroup() %>% 

  dplyr::select(-c(genus)) %>% 
  as.data.frame()

rownames(df) <- genus_matrix$genus
df_scaled <- as.data.frame(scale(df))

df_dist <- dist(df)

get_clust_tendency(df, n = 50)


# cluster method --------------------------------------------------------------

m <- c( "average", "single", "complete", "ward.D2")
names(m) <- c( "average", "single", "complete", "ward.D2")

#function to compute agglomerative coefficient
ac <- function(x) {
 a <- hclust(dist(df), method = x)
 coef.hclust(a)
}


#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)


# cluster numbers --------------------------------------------------------------

# index
n_clust2 <- NbClust(df, 
        min.nc = 2, max.nc = 15, 
        method = "ward.D2", 
        index =c("all"))


res_nbclust <- n_clust2$Best.nc
res_nbclust <- n_clust2$All.index

number_of_clusters <- 5

# hierachical clustering assessment --------------------------------------------
res.hc <- eclust(df, "hclust", k = number_of_clusters, graph = F)

final_data_hc <- cbind(df, cluster_hc = res.hc$cluster)


sil <- silhouette(res.hc$cluster, dist(df))
rownames(sil) <- rownames(df)

codes <- c("Australia-centered", 
                            "Holarctic-centered", 
                            "Africa-centered", 
                            "Tropical-Americas-centered", 
                           "Indo-Malesia-centered")

fviz_silhouette(sil, print.summary = T, label = T) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  scale_fill_manual(labels=codes, values= c( "1" = "#E73F74", 
                                            "2"= "#225522",
                                            "3" = "#F1CE63",  
                                            "4"= "#9467BD", 
                                            "5" = "#255e96")) +
  scale_color_manual(labels=codes, values= c( "1" = "#E73F74", 
                                             "2"= "#225522",
                                             "3" = "#F1CE63",  
                                             "4"= "#9467BD", 
                                             "5" = "#255e96"))


# Silhouette width of observations
sil_values <- res.hc$silinfo$widths[, 1:3]

sil_values %>% 
  group_by(cluster) %>% 
  summarise(
    mean = mean(sil_width), 
    median = median(sil_width), 
    sd = sd(sil_width), 
    coef = sd / mean * 100
            )
# Objects with negative silhouette
sil_values %>% 
  filter(sil_width < 0)



# bootstrapping clusters -------------------------------------------------------

dist_df <- dist(df)
df_pv <- t(df)



res.pv <- pvclust(df_pv, method.hclust = "ward.D2",
                  method.dist="euclidian", 
                  parallel = F, nboot = 100)

print(res.pv)


plot(res.pv, hang = -1, cex = 0.3, add = F, print.pv = "si", cex.pv=0.5, col.pv=c(si=4, au=2, bp=3, edge=NULL))
pvrect(res.pv, alpha=0.9)


groups_pv_clust <- cutree(res.pv$hclust, k = number_of_clusters)


### Final data 
final_data_boot <- cbind(final_data_hc, cluster_pvclust = groups_pv_clust) %>% 
  mutate(genus = rownames(.)) %>% 
  mutate(cluster = cluster_hc, 
         correspondence = cluster_hc == cluster_pvclust
  ) %>% 
 mutate(cluster = case_when(cluster == 1 ~ "Australia-centered", 
                            cluster == 2 ~ "Holarctic-centered", 
                            cluster == 3 ~ "Africa-centered", 
                            cluster == 4 ~ "Tropical-Americas-centered", 
                            cluster == 5 ~ "Indo-Malesia-centered"))



# plotting ---------------------------------------------------------------------

### Heatmap


col_names <- colnames(df)
clusters <- as.data.frame(final_data_boot$cluster)
colnames(clusters) <- "cluster"
rownames(clusters) <- rownames(df)
color_pal <- cols4all::c4a("friendly11", n = 9)



ann_colors = list(
  cluster =    c( "Australia-centered" = "#E73F74", 
                  "Holarctic-centered" = "#225522",
                  "Africa-centered" = "#F1CE63",  
                  "Tropical-Americas-centered"= "#9467BD", 
                  "Indo-Malesia-centered"= "#255e96") 
  
)

newnames <- lapply(
  rownames(clusters),
  function(x) bquote(italic(.(x))))


my_heatmap <- pheatmap(df, clustering_method = "ward.D2", cutree_rows = number_of_clusters, cluster_cols = F, 
                       cluster_rows = T, clustering_distance_rows = "euclidean", 
         color = viridis(n = 100, option = "D"), 
         annotation_row = clusters, 
         labels_row = as.expression(newnames), 
         #annotation_col = realm, 
         treeheight_row = 0, 
         treeheight_col = 0, 
         border_color = "black", 
         cellwidth = 20, 
         cellheight = 10, 
         annotation_colors = ann_colors, 
         show_rownames = TRUE)

color_pal <- cols4all::c4a("friendly11", n = 5)

### PCA hierachical 

df_pca_able <- final_data_boot
clusters <- df_pca_able[, "cluster"]


df_pca <- df_pca_able %>% 
  dplyr::select(-c(contains("cluster"), correspondence, genus))


hierachcial_pca <- fviz_pca_biplot(prcomp(df_pca, scale. = F),              # Visualize clusters in biplot
                        col.var = "black",
                        alpha.var = 0.6,
                        label = "all",
                        habillage = clusters, 
                        repel = TRUE,
                        addEllipses = TRUE,
                        ellipse.type = "convex",
                        #ellipse.level = 0.99, 
                        gradient.cols	 = "viridis", 
                        labelsize = 2, 
                        select.var = list(contrib=10)) +
  scale_color_manual(values =    c( "Australia-centered" = "#E73F74", 
                                                 "Holarctic-centered" = "#225522",
                                                 "Africa-centered" = "#F1CE63",  
                                                 "Tropical-Americas-centered"= "#9467BD", 
                                                 "Indo-Malesia-centered"= "#255e96") ) +
  scale_fill_manual(values =   c( "Australia-centered" = "#E73F74", 
                                                "Holarctic-centered" = "#225522",
                                                "Africa-centered" = "#F1CE63",  
                                                "Tropical-Americas-centered"= "#9467BD", 
                                                "Indo-Malesia-centered"= "#255e96") ) +
  scale_linetype_manual(values = c(1,1,1,1,1,1,1,1,0)) +
  theme_bw(base_size = 12)  



### Dendrogramm

phylo_tree <- as.phylo(res.hc)

group_colors <- c("1" = "#E73F74", 
                  "2" = "#225522",
                  "3" = "#F1CE63",  
                  "4" = "#9467BD", 
                  "5" = "#255e96")


cluster_assignments <- as.factor(cutree(res.hc, k = number_of_clusters))


tip_data <- data.frame(
  label = phylo_tree$tip.label,  
  cluster = cluster_assignments  
)



unrooted_dend <- ggtree(phylo_tree, layout = "unrooted") %<+% tip_data +
  geom_tippoint(aes(color = cluster),  size = 1, show.legend = F) +
  scale_color_manual(values =  group_colors) +

    geom_cladelabel(node=91, label = "", color="#E73F74", align = F) + 
    geom_cladelabel(node=92, label = "", color="#F1CE63", align = F) + 
    geom_cladelabel(node=89, label = "",  color="#255e96", align = F) + 
    geom_cladelabel(node=87, label = "",  color="#225522", align = F) + 
    geom_cladelabel(node=85, label = "",  color="#9467BD", align = F) + 
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(fill = NA)
  )


(dend_pca <- unrooted_dend + hierachcial_pca +  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold')))


