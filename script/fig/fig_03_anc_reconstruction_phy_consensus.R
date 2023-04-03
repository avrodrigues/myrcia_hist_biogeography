
# load packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(ape)
library(ggtree)
library(tidytree)
library(rnaturalearth)

l.func <- list.files(here("function"), full.names = T)
map(l.func, source)

source("../Shiny/map_clade_richness/R/functions.R")
source("../Shiny/map_clade_richness/R/gg_map_clade_richness.R")



# load data -------------------------------------------

# |- phy ----
phy.path <- here(
  "data", "phylogeny", "phy_cleaned", "000_phy_myrcia_cleaned_consensus.new"
)
myrcia_tree <- read.tree(phy.path)

# |- tipranges ----
geog.path <- here(
  "output", "biogeobears", "spp_area", "000_areas_myrcia_phy_consensus.data"
)
tipranges <- tipranges_from_phylip_tree(geog.path, myrcia_tree)

# |- model of recostruction ----
l_biogeo_mod <- readRDS(
  here("output",
       "biogeobears", 
       "list_results_models_biogeobears_phy_consenso.rds")
)

# best model 
l_biogeo_mod$table_AICc %>% 
  dplyr::arrange(AICc) %>% 
  mutate(delta_AICc =  AICc - min(AICc))

max_range_size = 2
# |- areas for nodes ----
areas_node <- 
  get_node_range_BioGeoBEARS(
    l_biogeo_mod$resDECj,
    geog.path,
    myrcia_tree,
    max_range_size
  )

# Figure -----------------------------------------------------


# |- prepare phy data for plot ----
fig_tree <- ggtree(myrcia_tree)

tree_df <- 
  fig_tree$data %>% 
  mutate(
    x = x-max(x),
    label = paste("M.", word(label, 2,2,"_")),
    area = c(tipranges, areas_node$area)
  ) 

tree_df <- tree_df %>% 
  mutate(
    parent_area = map_chr(node, function(x){ 
      if(x == 97) {
        return(tree_df %>% filter(node == 97) %>% pull(area))
      }
      parent(tree_df, x) %>% 
        pull(area)
    })
  )

# |- set range colors ----
unique(tree_df$area)

range_cols  <- 
  c( 
    "A" = "#37396C",
    "B" = "#58A0A9", 
    "C" = "#D3A838",
    "D" = "#362401",
    "E" = "#BB4455",
    
    "AB" = "#5A7C9F",
    #"ABC" = "#3B975D",
    #"AC" = "#4D7D38",
    "BC" =  "#3B975D",
    #"BCD" = "#75794B",
    #"BCE" = "#E4938F",
    
    "W" = "grey50"
  )

#show_colors(range_cols)

(g <-
  ggtree(
    tree_df, 
    layout = "rect",
    aes(x = x, color=parent_area), 
    size = 2
  )  +
  geom_tiplab(
    color = "#202015", 
    offset = 0.5, 
    fontface = 3
  ) +
  geom_label(
    aes(x = x, y = y, 
        label = area, #rep(" ", 191), 
        fill = area), 
    color = "#ffffff", 
    label.r = unit(0.01, "lines"), 
    label.padding = unit(0.2, "lines"), 
    size = 2.5, 
    fontface = 2
  ) +
  scale_fill_manual(name = "Areas", values = range_cols) +
  scale_color_manual(name = "Areas", values = range_cols) +
  scale_x_continuous(
    limits = c(-30, 5), 
    breaks = seq(-30, 5, 5), 
    labels = c(paste0(seq(30, 0, -5), "Ma"), "") 
  ) +
  theme_tree2() +
  theme(
    legend.position = "bottom", 
    legend.title = element_text(
      size = 18
    ),
    legend.text = element_text(
      size = 18
    ),
    axis.text = element_text(
      size = 18
    )
  ) 
)

ggsave(
  "output/fig/03_anc_rec_phy_consensus.png",
  g , 
  width = 10, 
  height = 15.8
) 


