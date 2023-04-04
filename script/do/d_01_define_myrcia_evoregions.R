
# load packages -----------------------------------------------------------
#devtools::install_github("GabrielNakamura/Herodotools")

library(tidyverse)
library(here)
library(ape)
library(furrr)
library(Herodotools)


# load data ---------------------------------------------------------------

# |- composition ----
comp_W <- read.csv(here("data", "W.csv"))
site_xy <-  read.csv(here("data", "W_xy.csv"))


# |- phylogeny ----
myrcia_tree_consensus <- read.tree(
  here(
    "data", 
    "phylogeny",
    "phy_cleaned",
    "000_phy_myrcia_cleaned_consensus.new")
)

# evoregions ----------------------------------------------------------


pcps_res <- PCPS::pcps(comp_W, cophenetic(myrcia_tree_consensus))
num_axis <- sum(pcps_res$values$Relative_eig > 0.05)
pcps_vectors <- pcps_res$vectors


n_max_evoreg <- find_max_nclust(
  pcps_vectors, 
  threshold = num_axis, 
  nperm = 1000,
  max.nclust = c(8, 10, 12, 14),
  subset = round(nrow(pcps_vectors)*0.1),
  confidence.level = 0.95
  )

#|- consensus phylogeny ----
l_evo_consensus <- evoregions(
  comp_W, 
  myrcia_tree_consensus, 
  max.n.clust = 8
)


save_dir_evo_consensus <- here("output", "evoregion", "phy_consensus")
save_file <- paste0("res_evoregion_phy_consensus.rds")

saveRDS(l_evo_consensus, here(save_dir_evo_consensus, save_file))
  
