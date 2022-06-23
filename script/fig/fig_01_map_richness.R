
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(rcartocolor)
library(sf)


# load local functions
l.func <- list.files(here("function"), full.names = T)
purrr::map(l.func, source)

# load map themes
source(here("script", "fig", "fig_00_map_themes.R"))

# loda data ---------------------------------------------------------------

# myrcia_dist was loaded from 'fig_00_map_themes.R'

# myrcia tree 
phy.path <- here(
  "data", "phylogeny", "phy_cleaned", "000_phy_myrcia_cleaned_consensus.new"
)
myrcia_tree <- read.tree(phy.path)
myrcia_tree_df <- generate_tree_df(myrcia_tree)


# map myrcia richness -----------------------------------------------------

rich_df <- 
  myrcia_dist %>% 
  count(x,y, name = "richness")

n_species <- length(unique(myrcia_dist$species))

map_richness <- 
ggplot() +
  geom_raster(data = rich_df, aes(x, y, fill = richness)) +
  labs(
    title = "*Myrcia* richness (n = 307)",
    fill = "Richness") +
  theme_map_continuous +
  theme(
    legend.background = element_rect(fill = NA)
  ) 
map_rich <- lemon::reposition_legend(map_richness, "bottom left", offset = 0.01)

mycia_clade_df <- select_clade(myrcia_tree_df, 97)

map_rich_tree <- map_clade_richness(
  mycia_clade_df, 
  myrcia_dist, 
  sf_countries, 
  color_pkg = "rcartocolors", 
  title = "Richness from species in phylogeny (n = 96)", 
  plot_theme = theme_map_continuous
) %>% 
  lemon::reposition_legend("bottom left", offset = 0.01)

chart_richness <- wrap_plots(
  list(map_rich, map_rich_tree),
  nrow = 1, ncol = 2
) +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

ggsave(
  here("output", "fig", "01_chart_richness.png"), 
  chart_richness, 
  width = 10,
  height = 5
)
