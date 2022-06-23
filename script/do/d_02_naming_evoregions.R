
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(rcartocolor)

source(here("function", "stand_evoregion_names.R"))

# load data ---------------------------------------------------------------

site_xy <- read.csv(here("data", "W_comp.csv"))

path_evo_consensus <- here(
  "output", "evoregion", "phy_consensus", "res_evoregion_few_occ_thr_site.rds"
)

path_evo_200 <- list.files(
  here("output", "evoregion", "phy_200"), 
  full.names = T
)

evoregions_df <- 
map_dfc(c(path_evo_consensus, path_evo_200), function(x){
  readRDS(x)$Cluster_Evoregions
}) %>%  
  set_names(
    c("phy_con", paste0("phy_", 1:200))
  ) %>% 
  as.data.frame()


# standardize evoregion names ---------------------------------------------

evo_ref <- c("C", "A", "E", "B", "D", "F", "G")

evo_df_std <- stand_evoregion_names(evoregions_df, evo_ref)  %>%  
  set_names(
    c("phy_con", paste0("phy_", 1:200))
  ) %>% 
  as.data.frame()

evo_df_std_xy <- bind_cols(site_xy, evo_df_std)

write.csv(
  evo_df_std_xy, 
  here("output", "evoregion", "evoregions_stand_names_df.csv"), 
  row.names = FALSE
)

evo_df_xy <- bind_cols(site_xy, evoregions_df)

not_five_evoregions <- 
apply(evoregions_df, 2, function(x){
  length(unique(x)) != 5
})  %>% 
  which() %>% 
  names()




#  |- check standardization -----------------------------------------------

# map theme

blue_gold_red_2 <- c(
  "#363870",
  "#589eab",
  "#d3a838",
  "#362401",
  "#cb3b2e",
  "#f9c0c2",
  "#346F3B"
)
names(blue_gold_red_2) <- LETTERS[1:7]

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
    title = element_text(color = greys[2]),
    axis.text = element_text(color = greys[2]), 
    axis.ticks = element_line(color = greys[3]), 
    panel.border = element_rect(color = greys[4], fill = NA), 
    axis.title = element_blank()
  )
)
coast <- rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)

map.limits <- list(
  x = c(-100, -30),
  y = c(-55, 25)
)

# map function
evo_plot <- function(evo_data, col_name){
  
  ggplot() +
    geom_raster(
      data = evo_data, 
      aes(x = x, y = y, fill = .data[[col_name]])
    ) +
    scale_fill_manual(
      name = "Evoregion",
      values = (blue_gold_red_2[]),
      na.value = "#C6C6C2"
    ) +
    geom_sf(data = coast) +
    coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
    ggtitle(col_name) +
    theme(
      legend.position = "bottom"
    ) +
    my_theme
}

map_evo_ref <- evo_plot(evo_df_std_xy, "phy_con")

# |-- generate map comparisons ----

for(i in seq_along(not_five_evoregions)){
  map_evo_compare <- evo_plot(evo_df_std_xy, not_five_evoregions[[i]])
  Sys.sleep(2)
  print(
    map_evo_ref + map_evo_compare
    )

}


# frequence of classification to the same evoregion -----------------------


freq_evo_reg <- site_freq_evoreg(evo_df_std[,-1])

freq_evo_reg$df.max.var.freq$evo.var.stand

evo_summary_df_xy <- bind_cols(site_xy, freq_evo_reg$df.max.var.freq)

map_evo_ref <- evo_plot(evo_df_std_xy, "phy_con")
map_evo_freq <- evo_plot(evo_summary_df_xy, "evo.names.max.freq")

map_evo_ref + map_evo_freq


sum(!evo_df_std_xy$phy_con == evo_summary_df_xy$evo.var.stand)

var_map <- 
ggplot() +
  geom_raster(
    data = evo_summary_df_xy, 
    aes(x = x, y = y, fill = .data[["evo.var.stand"]])
  ) +
  scale_fill_carto_c(
    palette = 7,
    #limits = c(0, 1)
  ) +
  #scale_fill_manual(
  #  name = "Evoregion",
  #  values = (blue_gold_red_2[]),
  #  na.value = "#C6C6C2"
  #) +
  geom_sf(data = coast) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  ggtitle("evo.var.stand") +
  theme(
    legend.position = "bottom"
  ) +
  my_theme

map_evo_ref + var_map
