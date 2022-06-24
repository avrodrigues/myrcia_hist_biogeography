
# load packages -----------------------------------------------------------

library(tidyverse)
library(vegan)
library(Herodotools)
library(here)


# load data ---------------------------------------------------------------

W <- read.csv(here("data", "W.csv"))
W_xy <- read.csv(here("data", "W_xy.csv"))

evo_df_std_xy <- read.csv(
  here("output", "evoregion", "evoregions_stand_names_df.csv")
)

regions <- readRDS(here("output/evoregion/phy_consensus/res_evoregion_few_occ_thr_site.rds"))

# Selecting only axis with more than 5% of explained variance from evoregion output
axis_sel <- which(regions$PCPS$prop_explainded >= regions$PCPS$tresh_dist)
PCPS_thresh <- regions$PCPS$vectors[, axis_sel] 

# distance matrix using 4 significant PCPS axis accordingly to the 5% threshold 
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# calculating affiliation values for each assemblage 
afi <- affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
                          groups = regions$Cluster_Evoregions) 

afilliation_data <- 
bind_cols(W_xy, afi) %>% 
  mutate(
    cat_afilliation = cut(
      afilliation, 5, 
      labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
      )
    )

map_afilliation <- 
ggplot(afilliation_data) +
  geom_raster(aes(x, y, fill = cat_afilliation)) +
  geom_sf(data = sf_coast, fill = NA) +
  coord_sf(xlim = evo_limits$x, ylim = evo_limits$y) +
  theme_evoregions +
  scale_fill_carto_d(
    name = "Affiliation",
    palette = "Geyser",
    direction = -1) +
  theme(
    legend.position = "bottom"
  )


ggsave(
  here("output", "fig", "s2_afilliation_evoregion.png"),
  map_afilliation, 
  width = 7,
  height = 7
)


