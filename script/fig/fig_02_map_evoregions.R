
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(sf)
library(terra)
library(rcartocolor)
library(patchwork)


source(here("function", "site_uncertainty.R"))
source(here("script/fig/fig_00_map_themes.R"))

# load data ---------------------------------------------------------------

evo_df_std_xy <- read.csv(
  here("output", "evoregion", "evoregions_stand_names_df.csv")
)


# prepare data ------------------------------------------------------------

evo_var <- site_uncertainty(evo_df_std_xy[, -c(1:3)])
evo_var_df_xy <- bind_cols(evo_df_std_xy[,1:2], evo_var = evo_var)

evo_limits <- list(
  x = range(evo_var_df_xy$x),
  y = range(evo_var_df_xy$y)
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


# evoregions map ----------------------------------------------------------

# |- consensus ----

map_evo_consensus <- 
ggplot() +
  geom_raster(
    data = evo_df_std_xy, 
    aes(x = x, y = y, fill = phy_con),
    #position = position_nudge(x = -0.25, y = 0.25)
  ) +
  scale_fill_manual(
    name = "Evoregion",
    values = blue_gold_red_2
    ) +
  geom_sf(
    data = sf_evoregion, 
    fill = NA, 
    color = bg, 
    size = 0.075) +
  geom_sf(data = sf_coast, fill = NA) +
  coord_sf(xlim = evo_limits$x, ylim = evo_limits$y) +
  theme(
    legend.position = "bottom"
  ) +
  theme_evoregions

#|- uncertainty ----

map_evo_uncertainty <- 
  ggplot() +
  geom_raster(
    data = evo_var_df_xy, 
    aes(x = x, y = y, fill = evo_var),
    #position = position_nudge(x = -0.25, y = 0.25)
  ) +
  scale_fill_carto_c(
    name = "Phylogenetic\nUncentainty",
    #colors = met.brewer("Greek", direction = -1)
    palette = "Sunset"
  ) +
  geom_sf(
    data = sf_evoregion, 
    fill = NA, 
    #color = bg, 
    size = 0.075) +
  geom_sf(data = sf_coast, fill = NA) +
  coord_sf(xlim = evo_limits$x, ylim = evo_limits$y) +
  theme(
    legend.position = "bottom"
  ) +
  theme_evoregions +
   guides(
     fill =  guide_colorbar(
       barheight = unit(5, "mm"), 
       label.hjust = 0,
       title.vjust = 1
       )
     )


map_evoregion <- map_evo_consensus + map_evo_uncertainty +
  plot_annotation(tag_levels = "A")

ggsave(
  here("output", "fig", "02_evoregion.png"),
  map_evoregion, 
  width = 9,
  height = 5
)
