
# load packages -----------------------------------------------------------
library(tidyverse)
library(here)
library(ape)
library(ggtree)
library(tidytree)
library(rnaturalearth)
library(Herodotools)

l.func <- list.files(here("function"), full.names = T)
map(l.func, source)

# source("../Shiny/map_clade_richness/R/functions.R")
# source("../Shiny/map_clade_richness/R/gg_map_clade_richness.R")



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

max_range_size = 3
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

    "BC" =  "#90c5a4",
    "BCD" = "#545c57",
    "BCE" =  "#88975D"
  )

tree_df$area <- factor(tree_df$area, levels = unique(tree_df$area))

(g <-
  ggtree(
    tree_df, 
    layout = "rect",
    aes(x = x, color= parent_area), 
    size = 2,
    show.legend = F
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
  scale_fill_manual(name = "Ranges", values = range_cols) +
  scale_color_manual(name = "Ranges", values = range_cols) +
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


# Supplement - Fig with pie charts in the nodes ----------------------------------------
theme_set(theme_bw())

range_cols_2  <- 
  c( 
    "A" = "#37396C",
    "B" = "#58A0A9", 
    "C" = "#D3A838",
    "D" = "#362401",
    "E" = "#BB4455",
    
    "AB" = "#5A7C9F",
    "AC" = "#3e9d57",
    
    "BC" =  "#90c5a4",
    "BCD" = "#545c57",
    "BD" = "#286c2e",
    
    "BCE" =  "#88975D",
    "BE" = "#674ab0",
    
    "other"  = "grey70"
  )

areas_node <- 
  get_node_range_BioGeoBEARS(
    l_biogeo_mod$resDECj,
    geog.path,
    myrcia_tree,
    max_range_size
  )

tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = geog.path)
trtable <- BioGeoBEARS::prt(myrcia_tree, printflag=FALSE)
areas <- BioGeoBEARS::getareas_from_tipranges_object(tipranges)
prob_state <- l_biogeo_mod$resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
states_list_0based <- cladoRcpp::rcpp_areas_list_to_states_list(areas = areas, 
                                                                maxareas = 3, 
                                                                include_null_range = TRUE)
ranges_list = NULL
for (i in 1:length(states_list_0based)) {
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

range_probabilities <- as.data.frame(prob_state)
row.names(range_probabilities) <- trtable$node
names(range_probabilities) <- ranges_list

prob_nodes <- range_probabilities |> 
  mutate(node = row_number(), .before = 1) |> 
  filter(!node %in% 1:96)

pie_chart_df <- 
apply(prob_nodes, 1, function(node){
  prob_vec <-  node[-1] |> sort(decreasing = T)
  
  n_most <- 3
  
  most_prob <- prob_vec[1:n_most]
  sum_other_prob <- sum(prob_vec[-c(1:n_most)])
  
  chart_df <- data.frame(
    ranges = c(names(most_prob), "other"),
    prob = c(as.numeric(most_prob), sum_other_prob),
    node = node[1]
  ) |> 
    arrange(desc(prob))
  
  chart_df
}) |> list_rbind() |> 
  pivot_wider(
    names_from = "ranges", 
    id_cols = "node", 
    values_from = prob,
    values_fill = 0
  )
  




pies <- nodepie(pie_chart_df, 2:ncol(pie_chart_df), outline.color = "white",  outline.size = .15)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = range_cols_2))




pie_tree <- ggtree(
  tree_df, 
  layout = "rect",
  aes(x = x), 
  color = "grey30",
  size = 2,
  show.legend = F
) +
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
  scale_fill_manual(name = "Ranges", values = range_cols_2) +
  scale_color_manual(name = "Ranges", values = range_cols_2) +
  geom_inset(pies, width = .055, height = .055) +
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


ggsave(
  "output/fig/S4_anc_rec_pie_chart.png",
  pie_tree, 
  width = 10, 
  height = 15.8
) 



# creates a legend plot

leg_df <- pie_chart_df[1,2:ncol(pie_chart_df)]


leg_df <- data.frame(
  area = factor(names(range_cols_2), levels = names(range_cols_2)), 
  value = runif(13, 6,8)
)

leg_plot <- ggplot(leg_df, aes(x = 1, y = value, fill = area)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(name = "Ranges", values = range_cols_2, 
                    labels = c(names(range_cols_2)[-13], "Other")) +
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

leg_only <- ggpubr::as_ggplot(ggpubr::get_legend(leg_plot))

ggsave(
  "output/fig/S4_legend_pie_chart.png",
  leg_only, 
  width = 10, 
  height = 2
) 
