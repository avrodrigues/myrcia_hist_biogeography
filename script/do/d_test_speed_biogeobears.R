# test the velocity in BiogeoBEARS

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

max_range_size = 4
num_cores_to_use = 8


# run with optmx ----------------------------------------------------------

start_optmix <- Sys.time()
l_biogeo_mod <- model_biogeobears_test(
      phy.path[[1]],
      geog.path[[1]],
      max_range_size,
      num_cores_to_use, 
      setup_optimx = TRUE
    )
end_optmix <- Sys.time()

time_spend_optimx <- end_optmix - start_optmix


# run with GenSA ----------------------------------------------------------


start_gensa <- Sys.time()
l_biogeo_mod <- model_biogeobears_test(
  phy.path[[1]],
  geog.path[[1]],
  max_range_size,
  num_cores_to_use, 
  setup_optimx = "GenSA"
)
end_gensa <- Sys.time()

time_spend_gensa <- end_gensa - start_gensa
