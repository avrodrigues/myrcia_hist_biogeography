
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(ape)
library(ggtree)
library(tidytree)
library(rnaturalearth)
library(rcartocolor)
library(lemon)
library(ggtext)
library(glue)

l.func <- list.files(here("function"), full.names = T)
map(l.func, source)


# load data ---------------------------------------------------------------

sf_coast <- ne_coastline(scale = 50, returnclass = "sf")
sf_countries <- ne_countries(scale = 50, returnclass = "sf")
sf_evoregion <- readRDS(here("output", "evoregion", "sf_evoregion.rds"))

# |- occ data ----
myrcia_occ <- readRDS("data/distribution/list_myrcia_binary_df_05_degree.rds")[[3]]

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

# |- areas for nodes ----
areas_node <- 
  get_node_range_BioGeoBEARS(
    l_biogeo_mod$resBAYAREALIKEj,
    geog.path,
    myrcia_tree,
    max_range_size
  )

# Figure -----------------------------------------------------

# |- set range colors ----

range_cols  <- 
  c( 
    "A" = "#37396C",
    "AB" = "#496B8E",
    "B" = "#58A0A9", 
    
    "ABCE" = "#72934A",
    "BCE" = "#9FC077",
    "CE" = "#FFB545",
    "C" = "#D3A838",
    
    "AE" = "#8F348A",
    "ABE" = "#B459AF",
    "BE" =  "#E589B6",
    
    "ABDE" = "#734358",
    "BDE" =  "#855573",
    
    "DE" = "#8F2400",
    "E" = "#CD3B2E"
  )

carto_names <- cartocolors %>% 
  filter(Type == "quantitative") %>% 
  pull(Name)

greys <- c(
  "#040400",
  "#1F1F1B",
  "#3B3B37",
  "#575753",
  "#73736F",
  "#8E8E8A",
  "#AAAAA6",
  "#C6C6C2",
  "#E2E2DE",
  "#FEFEFA"
)

my_theme <- list(
  theme(
    panel.background = element_rect(fill = "#FAF8F4"), 
    panel.grid = element_blank(), 
    text = element_text(color = greys[1]), 
    plot.title = element_markdown(color = greys[2]),
    axis.text = element_text(color = greys[2]), 
    axis.ticks = element_line(color = greys[3]), 
    panel.border = element_rect(color = greys[4], fill = NA), 
    axis.title = element_blank(), 
    legend.background = element_rect(fill = NA), 
    legend.text = element_text(size = 7), 
    legend.title = element_text(face = "bold", size = 9, margin = margin(b = 2)), 
  )
)

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

# exploring clade richness ------------------------------------------------

myrcia_tree_df <- generate_tree_df(myrcia_tree) %>% 
  left_join(tree_df %>% select(-label)) 

interactive_phy(myrcia_tree_df)

clade_names_df <- tribble(
  ~node, ~clade_name, 
  104, "Gomidesia",
  127, "Aguava",
  136, "Eugeniopsis",
  141, "Tomentosae",
  #102, "Gomideia + Agauva + Eugeniopsis + Tomentosae",
  142, "Clade 10",
  144, "Calyptramthes",
  153, "Sympodiomyrcia",
  160, "Reticulosae",
  164, "Myrcia",
  171, "Aulomyrcia"
)


list_clade_maps <- map(
  1:nrow(clade_names_df), 
  function(i){
    sel_clade_df <- select_clade(
      myrcia_tree_df, 
      clade_names_df$node[i]
    ) 
    
    n_tips <- 
      sel_clade_df %>% 
      filter(node_selected, isTip) %>% 
      nrow()
    
    clade <- clade_names_df$clade_name[i]
      
    title <- glue("{clade} (n = {n_tips})")
    
    plot_map <- map_clade_richness(
      sel_clade_df, 
      myrcia_occ, 
      sf_countries, 
      color_pkg = "rcartocolors", 
      title = title, 
      plot_theme = my_theme
    )
   
    reposition_legend(plot_map, "bottom left", offset = 0.01, plot = F)
  }
)

chart_clade_richness <-
  patchwork::wrap_plots(list_clade_maps, 
                      nrow = 5, ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

ggsave(
  here("output", "fig", "s1_clade_richness.png"), 
  chart_clade_richness, 
  width = 10,
  height = 20
)
