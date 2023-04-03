path_evo_consensus <- here(
  "output", "evoregion", "phy_consensus", "res_evoregion_few_occ_thr_site.rds"
)

site_xy <- read.csv(here("data", "W_xy.csv"))

evo <- readRDS(path_evo_consensus)

sum(evo$PCPS$prop_explainded > evo$PCPS$tresh_dist)

pcps_evo <- evo$PCPS$vectors[,1:3] |> vegan::decostand(method = "range")

col_vec <- 
apply(pcps_evo, 1, function(row){
  rgb(row[3], row[2], row[1], maxColorValue = 1)
})

unique_col <- unique(col_vec)
names(unique_col) <- unique_col

phylobeta_df <- site_xy |> mutate(col_vec)

ggplot(phylobeta_df, aes(x, y, fill = col_vec)) +
  geom_raster(show.legend = F) +
  scale_fill_manual(values = unique_col) +
  theme_bw()


## altitude --------------------------------------------

alt <- terra::rast(here("data","Altitude_neotropico.tif"))
alt_df <- as.data.frame(alt, xy = T)

library(rcartocolor)

alt_df_evo <- terra::extract(alt, as.matrix(site_xy), xy=TRUE) %>% 
  mutate(
    evoregion = evo$Cluster_Evoregions, 
    evo_names = c("C", "A", "E", "B", "D")[evoregion])

ggplot(alt_df_evo, aes(x = evo_names, y = log10(Altitude_neotropico))) +
  geom_boxplot()

ggplot(alt_df_evo, aes(x, y, fill = evo_names)) +
  geom_raster(show.legend = T) +
  theme_bw()

ggplot(alt_df_evo, aes(x, y, fill = sqrt(Altitude_neotropico))) +
  geom_raster(show.legend = T) +
  # scale_fill_gradientn(
  #   colors = terrain.colors(n = length(unique(sqrt(alt_df$Altitude_neotropico))))) +
  scale_fill_distiller(palette = "Spectral") +
  coord_equal() +
   
  theme_bw() + 
  theme( 
    legend.position = "bottom")
