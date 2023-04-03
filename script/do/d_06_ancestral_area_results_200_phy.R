
# load packages -----------------------------------------------------------

library(here)
library(tidyverse)

l.func <- list.files(here("function"), full.names = T)
walk(l.func, source)

# load data ---------------------------------------------------------------

l.bgb.files <- list.files(
  here("output", "biogeobears", "results_biogeobears_200_phy"), 
  full.names = T
  )

l.res <- map(l.bgb.files, read_rds)


phy.path <- here("data", "phylogeny", "phy_cleaned") |> 
  list.files(pattern = "posterior", full.names = T)

geog.path <- here("output", "biogeobears", "spp_area") |> 
  list.files(pattern = "posterior", full.names = T)

# best model --------------------------------------------------------------

best.mod.df <- map_df(seq_along(l.res), function(i){
  
  l.res[[i]]$table_AICc %>% 
    dplyr::arrange(AICc) %>% 
    mutate(delta_AICc =  AICc - min(AICc)) %>% 
    slice(1) %>% 
    mutate(
      best_model = row.names(.),
      phy = i
    )
  
})

# results from best model -------------------------------------------------

areas_node <- 
  get_node_range_BioGeoBEARS(
    l_biogeo_mod$resBAYAREALIKEj,
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
  l_biogeo_mod$resBAYAREALIKEj
)
dev.off()


get_node_range_BioGeoBEARS_v2(l_biogeo_mod$resBAYAREALIKEj)


geogfn <- l_biogeo_mod$resBAYAREALIKEj$inputs$geogfn
max_range_size <- l_biogeo_mod$resBAYAREALIKEj$inputs$max_range_size
include_null_range <- l_biogeo_mod$resBAYAREALIKEj$inputs$include_null_range

relprobs_matrix <- l_biogeo_mod$resBAYAREALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node


tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
areas = getareas_from_tipranges_object(tipranges)

statenames = areas_list_to_states_list_new(areas, 
                                           maxareas = max_range_size, 
                                           include_null_range = include_null_range, 
                                           split_ABC = FALSE)


MLstates = get_ML_states_from_relprobs(relprobs_matrix, 
                                       statenames, 
                                       returnwhat = "states")

n_sp <- nrow(tipranges@df)

png(
  filename = "output/fig/plot_ancstates_rascunho.png",
  height = 1400,
  width = 700, 
  units = "px"
  
)
plot(myrcia_tree, cex = 0.5)
nodelabels(MLstates[-c(1:n_sp)], cex = 0.5)
tiplabels(MLstates[c(1:n_sp)], cex = 0.5)
axisPhylo()
dev.off()





