# load packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(ape)

l.func <- list.files(here("function"), full.names = T)
map(l.func, source)


# ancentral area reconstruction -------------------------------------------

phy.path <- here("data", "phylogeny", "phy_cleaned") |> 
  list.files(pattern = "posterior", full.names = T)

geog.path <- here("output", "biogeobears", "spp_area") |> 
  list.files(pattern = "posterior", full.names = T)


# save path
one200 <- as.character(1:200)

one200_vec <- ifelse(
  nchar(one200) == 1, 
  paste0("00", one200), 
  ifelse(
    nchar(one200) == 2, 
    paste0("0", one200),
    one200)
)

dir_save <- here("output", "biogeobears", "results_biogeobears_200_phy")



max_range_size = 4
num_cores_to_use = 3

library(furrr)

plan(multisession, workers = 3)

future_map(1:200, function(i){
  l_biogeo_mod <- model_biogeobears(
    phy.path[[i]],
    geog.path[[i]],
    max_range_size,
    num_cores_to_use
  )
  
  path_save <- here(
    dir_save, paste0(one200_vec[i],"_phy_list_results_biogeobears.rds") 
  )
  
  saveRDS(l_biogeo_mod, file = path_save)
})


