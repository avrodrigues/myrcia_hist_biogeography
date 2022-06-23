
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(sf)
library(terra)
library(rnaturalearth)
library(MetBrewer)

source(here("function", "site_uncertainty.R"))

# load data ---------------------------------------------------------------

evo_df_std_xy <- read.csv(
  here("output", "evoregion", "evoregions_stand_names_df.csv")
)


# prepare data ------------------------------------------------------------

evo_var <- site_uncertainty(evo_df_std_xy[, -c(1:3)])
evo_var_df_xy <- bind_cols(evo_df_std_xy[,1:2], evo_var = evo_var)
  

# map theming -------------------------------------------------------------

# |- color option ----
blue_gold_red <- c(
  "#303260",
  "#88beca",
  "#e4b434",
  "#7b5a28",
  "#e19296",
  "#c62f22"
  
)
blue_gold_red_2 <- c(
  "#363870",
  "#589eab",
  "#d3a838",
  "#362401",
  "#cb3b2e",
  "#f9c0c2"
)

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

bg <- "#FAF8F4"

# |- ggplot theme ----
my_theme <- list(
  theme(
    panel.background = element_rect(fill = bg), 
    panel.grid = element_blank(), 
    text = element_text(color = greys[1]), 
    title = element_text(color = greys[2]),
    axis.text = element_text(color = greys[2]), 
    axis.ticks = element_line(color = greys[3]), 
    panel.border = element_rect(color = greys[4], fill = NA), 
    axis.title = element_blank()
  )
)

# |- continent and limits ----
coast <- rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)

map.limits <- list(
  x = c(-100, -30),
  y = c(-55, 25)
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

saveRDS(sf_evoregion, here("output", "evoregion", "sf_evoregion.rds"))

coast <- st_transform(coast, st_crs(sf_evoregion))


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
  geom_sf(data = coast, fill = NA) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme(
    legend.position = "bottom"
  ) +
  my_theme

#|- uncertainty ----

map_evo_uncertainty <- 
  ggplot() +
  geom_raster(
    data = evo_var_df_xy, 
    aes(x = x, y = y, fill = evo_var),
    #position = position_nudge(x = -0.25, y = 0.25)
  ) +
  scale_fill_gradientn(
    name = "Phylogenetic\nUncentainty",
    #colors = met.brewer("Greek", direction = -1)
    colors = met.brewer("OKeeffe2")
  ) +
  geom_sf(
    data = sf_evoregion, 
    fill = NA, 
    #color = bg, 
    size = 0.075) +
  geom_sf(data = coast, fill = NA) +
  coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
  theme(
    legend.position = "bottom"
  ) +
  my_theme +
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
  here("output", "fig", "01_evoregion_v2.png"),
  map_evoregion, 
  width = 9,
  height = 6.5
)
