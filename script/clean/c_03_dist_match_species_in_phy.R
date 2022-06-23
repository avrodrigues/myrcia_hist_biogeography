
# load packages -----------------------------------------------------------

library(ape)
library(tidyverse)


# load data ---------------------------------------------------------------

l_myrcia_binary_df <- read_rds(
  here("data", "distribution", "list_myrcia_binary_df_05_degree.rds")
)

myrcia_tree_consensus <- read.tree(
  here(
    "data", 
    "phylogeny",
    "phy_cleaned",
    "000_phy_myrcia_cleaned_consensus.new")
)


# organize data -----------------------------------------------------------

phy_species <- myrcia_tree_consensus$tip.label

l_myrcia_dist_phy <- map(
  l_myrcia_binary_df, function(x){
    
    x %>% 
      filter(species %in% phy_species) 
    
})


l_myrcia_dist_phy_wide <- map(
  l_myrcia_dist_phy, function(x){
    
    x %>% 
      pivot_wider(
        names_from = species, 
        values_from = presence, 
        values_fill = 0
      )
  }
)


l_myrcia_phy_site_xy <- map(
  l_myrcia_dist_phy_wide, function(x) x[,1:2]
  )
l_myrcia_phy_comp <-  map(
  l_myrcia_dist_phy_wide, function(x) x[,-c(1:2)]
)

# save data ---------------------------------------------------------------


write_rds(
  l_myrcia_dist_phy, 
  here("data", "distribution", "list_myrcia_dist_phy.rds")
)

write_rds(
  l_myrcia_phy_site_xy, 
  here("data", "distribution", "list_myrcia_phy_site_xy.rds")
)

write_rds(
 l_myrcia_phy_comp, 
  here("data", "distribution", "list_myrcia_phy_comp.rds")
)

