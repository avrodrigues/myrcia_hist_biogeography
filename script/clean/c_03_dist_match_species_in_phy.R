
# load packages -----------------------------------------------------------

library(ape)
library(tidyverse)


# load data ---------------------------------------------------------------

myrcia_binary_df <- read.csv(
  here("data", "distribution", "myrcia_binary_df_05_degree.csv")
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

myrcia_phy_long_format <- 
  myrcia_binary_df %>% 
      filter(species %in% phy_species) 
    
myrcia_phy_wide <- myrcia_phy_long_format %>% 
      pivot_wider(
        names_from = species, 
        values_from = presence, 
        values_fill = 0
      )

myrcia_phy_site_xy <- myrcia_phy_wide[,1:2]

myrcia_phy_comp <- myrcia_phy_wide[,-c(1:2)]

# save data ---------------------------------------------------------------


write.csv(
  myrcia_phy_long_format, 
  here("data", "distribution", "myrcia_phy_long_format.csv"), 
  row.names = F
)

write.csv(
  myrcia_phy_site_xy, 
  here("data", "distribution", "myrcia_phy_site_xy.csv"), 
  row.names = F
)

write.csv(
  myrcia_phy_comp, 
  here("data", "distribution", "myrcia_phy_comp.csv"), 
  row.names = F
)

