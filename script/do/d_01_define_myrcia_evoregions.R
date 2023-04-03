
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(ape)
library(furrr)
library(Herodotools)

source(here("function", "evoregions2.R"))

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

# l_phy_files <- list.files(
#   here("data", "phylogeny","phy_cleaned"), full.names = T
#   )[-1]

# myrcia_tree_200 <- l_phy_files %>% map(read.tree)


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
  
# #|- 200 phylogenies (posterior sampling) ----
# 
# 
# 
# cent <- rep(0:2, each = 100)[1:201]
# dec <- rep(rep(0:9, each = 10), 3)[1:201]
# unid <- rep(0:9, 21)[1:201]
# 
# um200 <- paste0(cent, dec, unid)[-1]
# 
# save_dir_evo_200 <- here("output", "evoregion", "phy_200")
# 
# plan(multisession, workers = 2)
# 
# furrr::future_map(
#   seq_along(myrcia_tree_200), 
#   function(i){
#     
#     evo <-  evoregions2(
#       comp_W, 
#       method.clust = "ward",
#       myrcia_tree_200[[i]], 
#       max.n.clust = 10
#     )
#     
#     file_name <- paste0("res_evoregion_phy_", um200[i], ".rds")
#     
#     saveRDS(evo, here(save_dir_evo_200, file_name))
#   }
# )
# 
