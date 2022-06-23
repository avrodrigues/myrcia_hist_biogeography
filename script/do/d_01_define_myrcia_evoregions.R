
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(ape)

source(here("function", "evoregions.R"))

# load data ---------------------------------------------------------------

# |- distribution ----

l_myrcia_phy_site_xy <- read_rds( 
  here("data", "distribution", "list_myrcia_phy_site_xy.rds")
)

l_myrcia_phy_comp <- read_rds(
  here("data", "distribution", "list_myrcia_phy_comp.rds")
)

# keep only sites with more than 2 species
l_keep_site_3sp <- map(l_myrcia_phy_comp, function(x){
  rowSums(x) >= 3
})
  

l_m_comp <- map(1:3, function(i){
  l_myrcia_phy_comp[[i]] %>% filter(l_keep_site_3sp[[i]])
})

write.csv(l_m_comp[[3]], here("data", "W.csv"), row.names = F)

l_m_site_xy <- map(1:3, function(i){
  l_myrcia_phy_site_xy[[i]] %>% filter(l_keep_site_3sp[[i]])
})

write.csv(l_m_site_xy[[3]], here("data", "W_comp.csv"), row.names = F)

# |- phylogeny ----
myrcia_tree_consensus <- read.tree(
  here(
    "data", 
    "phylogeny",
    "phy_cleaned",
    "000_phy_myrcia_cleaned_consensus.new")
)

l_phy_files <- list.files(
  here("data", "phylogeny","phy_cleaned"), full.names = T
  )[-1]


myrcia_tree_200 <- l_phy_files %>% map(read.tree)


# evoregions ----------------------------------------------------------

#|- consensus phylogeny ----
l_evo_consensus <- map(
  1:3, function(i){
    evo <- evoregions(
      l_m_comp[[i]], 
      myrcia_tree_consensus, 
      max.n.clust = 10
    )
    
    cat(paste("evoregion comp mtx", i), fill = T)
    return(evo)
  }
)
  
names(l_evo_consensus) <- names(l_myrcia_phy_site_xy)

save_dir_evo_consensus <- here("output", "evoregion", "phy_consensus")

for(i in seq_along(l_evo_consensus)) {
  save_file <- paste0("res_evoregion_", names(l_evo_consensus)[i], ".rds")
  
  saveRDS(l_evo_consensus[[i]], here(save_dir_evo_consensus, save_file))
  
}

#|- 200 phylogenies (posterior sampling) ----
evo_200 <- map(
  seq_along(myrcia_tree_200), 
  function(i){
    evo <-  evoregions(
      l_m_comp[[3]], 
      myrcia_tree_200[[i]], 
      max.n.clust = 10
    )
    
    cat(paste("evoregion phy", i), fill = T)
    return(evo)
  }
)

