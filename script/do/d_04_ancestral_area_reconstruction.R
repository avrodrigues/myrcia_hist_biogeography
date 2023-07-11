
# load packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(ape)

l.func <- list.files(here("function"), full.names = T)
map(l.func, source)


# ancentral area reconstruction -------------------------------------------

phy.path <- here(
  "data", "phylogeny", "phy_cleaned", "000_phy_myrcia_cleaned_consensus.new"
  )
myrcia_tree <- read.tree(phy.path)

geog.path <- here(
  "output", "biogeobears", "spp_area", "000_areas_myrcia_phy_consensus.data"
  )

max_range_size = 3
num_cores_to_use = 10

l_biogeo_mod <- model_biogeobears(
    phy.path,
    geog.path,
    max_range_size,
    num_cores_to_use
)

saveRDS(l_biogeo_mod, "output/biogeobears/list_results_models_biogeobears_phy_consenso.rds")

# best model 
l_biogeo_mod$table_AICc %>% 
  dplyr::arrange(AICc) %>% 
  mutate(delta_AICc =  AICc - min(AICc))


# results from best model -------------------------------------------------

areas_node <- 
get_node_range_BioGeoBEARS(
  l_biogeo_mod$resDECj,
  geog.path,
  myrcia_tree,
  max_range_size
)

tipranges <- tipranges_from_phylip_tree(geog.path, myrcia_tree)

tipranges_mtx <- read_phylip_file(geog.path)
tip_mtx_order <- match(myrcia_tree$tip.label, rownames(tipranges_mtx))

tipranges_mtx <- tipranges_mtx[tip_mtx_order, ]



# plot ancestral reconstruction -------------------------------------------

#|- draft ----
plot(myrcia_tree, cex = 0.5)
nodelabels(areas_node$area, cex = 0.5)

png(
  filename = "output/fig/plot_ancstates_rascunho_BGB.png",
  height = 1400,
  width = 700, 
  units = "px"
  
)

plot_BioGeoBEARS_results(
  l_biogeo_mod$resDECj
)

dev.off()


