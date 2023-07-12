
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(sf)
library(terra)
library(rcartocolor)
library(patchwork)
library(ggtree)

# 
# source(here("function", "site_uncertainty.R"))
source(here("script/fig/fig_00_map_themes.R"))

# load data ---------------------------------------------------------------

evo_df_std_xy <- read.csv(
  here("output", "evoregion", "evoregions_stand_names_df.csv")
)

afi <- readRDS(
  here("output", "evoregion", "phy_consensus", "evoregion_afiliation.rds")
)

path_evo_consensus <- here(
  "output", "evoregion", "phy_consensus", "res_evoregion_phy_consensus.rds"
)
evo_consensus <- readRDS(path_evo_consensus)


# prepare data ------------------------------------------------------------

# evo_var <- site_uncertainty(evo_df_std_xy[, -c(1:3)])
evo_afi_df_xy <- bind_cols(
  evo_df_std_xy[,1:2], evo_afiliation = afi[,"afilliation"]
  ) %>% 
  mutate(
    cat_afilliation = cut(
      evo_afiliation, 5, 
      labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
    )
  )

evo_limits <- list(
  x = range(evo_afi_df_xy$x),
  y = range(evo_afi_df_xy$y)
)

# |- evoregion bondaries ----
evo_df <- data.frame(
  evo_df_std_xy[,1:2],
  phy_con = as.numeric(as.factor(evo_df_std_xy[,3]))
)

r_evoregion <- rast(evo_df)
sf_evoregion <- as.polygons(r_evoregion) %>% 
  st_as_sf()

st_crs(sf_evoregion) <- 4326

#saveRDS(sf_evoregion, here("output", "evoregion", "sf_evoregion.rds"))

# |- relationships among evoregions ----

# |- defining the reference names for evoregions
pcps_vsc_thr <- evo_consensus$PCPS$vectors[,1:3]
k_grps_df <- data.frame(evo_grp = evo_df_std_xy[,3], pcps_vsc_thr)


dist_mtx <- k_grps_df %>% 
  group_by(evo_grp) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  select(-evo_grp) %>% 
  dist() %>% 
  as.matrix()

dimnames(dist_mtx) <- list(LETTERS[1:5], LETTERS[1:5])

dendro_evo <- dist_mtx %>% 
  as.dist %>% 
  hclust("ward.D2") %>% 
  as.dendrogram() 

# evoregions map ----------------------------------------------------------

# |- consensus ----

(map_evo_consensus <- 
ggplot() +
  #geom_sf(data = sf_pol_coast , fill = greys[5]) +
  geom_raster(
    data = evo_df_std_xy, 
    aes(x = x, y = y, fill = phy_con),
    #position = position_nudge(x = -0.25, y = 0.25)
  ) +
  scale_fill_manual(
    name = "Evoregion",
    values = colors_evo
    ) +
  geom_sf(
    data = sf_evoregion, 
    fill = NA, 
    color = bg, 
    size = 0.13) +
  geom_sf(data = sf_coast, fill = NA) +
  coord_sf(xlim = evo_limits$x, ylim = evo_limits$y) +
  theme(
    legend.position = "bottom"
  ) +
  theme_evoregions
)

# |-- relationship among evoregions ----

dendro_df <- ggtree(dendro_evo)

dendro_df <- 
  dendro_df$data |> 
  mutate(
    x = ifelse(isTip, 0, x)
  )

dendro_plot <-  ggtree(dendro_df) +
  layout_dendrogram() +
  geom_label(
    aes(x = x, y = y, 
        label = label, 
        fill = label), 
    color = "#ffffff", 
    label.padding = unit(0.35, "lines"), 
    size = 3.5, 
    fontface = 2, 
    show.legend = F
  ) +
  scale_fill_manual(values = blue_gold_red_2) +
  theme_dendrogram() +
  theme(
    plot.background = element_rect(fill = NA, color = greys[5]),
    panel.background = element_rect(fill = NA, colour = NA))

dendro_plot <- ggtree::flip(dendro_plot, 5, 7)

# |-- Evoregion chart

chart_evo_dendro <- map_evo_consensus +
  inset_element(dendro_plot,
                clip = F,
                0.025, 0.025, 0.3, 0.35) +
  plot_annotation(tag_levels = "A")

ggsave(
  "output/fig/02_evoregion_chart.png", 
  chart_evo_dendro,
  height = 8,
  width = 8
)

#|- affiliation ----


(map_evo_affiliation <- 
  ggplot(evo_afi_df_xy) +
  geom_raster(aes(x, y, fill = cat_afilliation)) +
  geom_sf(data = sf_coast, fill = NA) +
  scale_fill_carto_d(
    name = "Affiliation",
    palette = "Geyser",
    direction = -1) +
  geom_sf(data = sf_coast, fill = NA) +
  coord_sf(xlim = evo_limits$x, ylim = evo_limits$y) +
  theme(
    legend.position = "bottom"
  ) +
  theme_evoregions 
)

ggsave(
  here("output", "fig", "s1_site_affiliation.png"),
  map_evo_affiliation, 
  width = 8,
  height = 8
)
