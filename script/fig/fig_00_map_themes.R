

# load packages -----------------------------------------------------------

library(rcartocolor)
library(ggtext)
library(ggplot2)


# load data ---------------------------------------------------------------

myrcia_dist <- read.csv(
  here("data", "distribution", "myrcia_binary_df_05_degree.csv")
)


# |- continent and limits ----
sf_coast <- rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)
sf_countries <- rnaturalearth::ne_countries(returnclass = "sf", scale = 50)

map.limits <- list(
  x = range(myrcia_dist$x),
  y = range(myrcia_dist$y)
)

# |- colors ----

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

colors_evo <-  c(
  "#363870",
  "#589eab",
  "#d3a838",
  "#362401",
  "#BB4455"
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


# |- ggplot themes ----

theme_map_continuous <- list(
  scale_fill_carto_c(palette = "Emrld", name = "Richness"),
  geom_sf(data = sf_countries, fill = NA),
  coord_sf(xlim = map.limits$x, ylim = map.limits$y),
  theme(
    panel.background = element_rect(fill = bg), 
    panel.grid = element_blank(), 
    text = element_text(color = greys[1]), 
    plot.title = element_markdown(color = greys[2]),
    axis.text = element_text(color = greys[2]), 
    axis.ticks = element_line(color = greys[3]), 
    panel.border = element_rect(color = greys[4], fill = NA), 
    axis.title = element_blank(), 
    legend.background = element_rect(fill = NA), 
    legend.text = element_text(size = 7), 
    legend.title = element_text(face = "bold", size = 9, margin = margin(b = 2))
  )
)


theme_evoregions <- list(
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
