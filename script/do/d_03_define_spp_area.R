
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)

source(here("function", "tipranges_to_BioGeoBEARS.R"))

# load data ---------------------------------------------------------------

# evoregion classification
evo_df <- read.csv(here("output", "evoregion", "evoregions_stand_names_df.csv"))
site_xy <- evo_df[,1:2]
evo_mtx <- evo_df[, -c(1:2)]

# species composition
myrcia_comp <- read.csv(here("data", "W.csv"))


# define species areas ----------------------------------------------------


  area_names <- unique(evo_mtx)
  area_names <- area_names[order(area_names)]
  
  phy_area_mtx <- matrix(
    0, 
    nrow = ncol(myrcia_comp), 
    ncol = length(area_names), 
    dimnames = list(
      names(myrcia_comp), 
      area_names)
    )
  
  for(sp in 1:ncol(myrcia_comp)){
    
    evo_spp <- evo_mtx[myrcia_comp[[sp]] == 1] 
    
    areas <- table(evo_spp)
    areas_prop <- areas/sum(areas)
    areas_dist <- names(areas_prop)[which(areas_prop > 0.45)]
    area_pos <- match(areas_dist, area_names)
    
    phy_area_mtx[sp, area_pos] <- 1
  }

  

  phy_area_mtx <-
  phy_area_mtx %>%
    as.data.frame %>%
    mutate(W = ifelse(rowSums(phy_area_mtx) == 0, 1, 0)) %>%
    as.matrix()


# save areas to biogeobears -----------------------------------------------

dir_save <- here("output", "biogeobears", "spp_area")

area_mtx <- phy_area_mtx# l_spp_area[[i]]
filename <- paste0("000_areas_myrcia_phy_consensus.data")

file_save <- here(dir_save, filename)
  
tipranges_to_BioGeoBEARS(
    area_mtx, 
   filename = file_save, 
   areanames = colnames(area_mtx)
  )



