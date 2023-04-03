
# load packages -----------------------------------------------------------

library(here)
library(dplyr)

# load data ---------------------------------------------------------------

myrcia_phy_site_xy <- 
  read.csv(here("data", "distribution", "myrcia_phy_site_xy.csv"))

myrcia_phy_comp <- 
  read.csv(here("data", "distribution", "myrcia_phy_comp.csv"))


# filter and save data ----------------------------------------------------

# keep only sites with at least 3 species
keep_site_3sp <- rowSums(myrcia_phy_comp) >= 3

#  composition 
m_comp <- myrcia_phy_comp %>% filter(keep_site_3sp)
write.csv(m_comp, here("data", "W.csv"), row.names = F)

# geographycal coordinates

site_xy <- myrcia_phy_site_xy %>% filter(keep_site_3sp)
write.csv(site_xy, here("data", "W_xy.csv"), row.names = F)
