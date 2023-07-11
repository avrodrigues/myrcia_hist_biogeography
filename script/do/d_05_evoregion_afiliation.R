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

regions <- readRDS(
  here("output", "evoregion", "phy_consensus", "res_evoregion_phy_consensus.rds")
  )

# Selecting only axis with more than 5% of explained variance from evoregion output
axis_sel <- which(regions$PCPS$prop_explainded >= regions$PCPS$tresh_dist)
PCPS_thresh <- regions$PCPS$vectors[, axis_sel] 

# distance matrix using 3 significant PCPS axis accordingly to the 5% threshold 
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean")

# calculating affiliation values for each assemblage 
afi <- calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
                          groups = regions$cluster_evoregions) 

saveRDS(
  afi, 
  here("output", "evoregion", "phy_consensus", "evoregion_afiliation.rds")
)
