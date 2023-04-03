library(biscale)
library(cowplot)
library(rcartocolor)
library(viridis)
library(patchwork)


# load map themes
source(here("script", "fig", "fig_00_map_themes.R"))

blue_grad <- c("#D3D3D3", "#98B8C0", "#5B9CAD")
pink_grad <- c("#D3D3D3", "#C59595", "#B65252")

# maps -------------------------------------------------

sites_db <- readRDS("sites_db_metrics.rds")

map_insitu <- 
sites_db %>% 
  ggplot() + 
  geom_raster(aes(x = x, y = y, fill = prop_insitu)) +
  # scale_fill_viridis_c(
  #   name = "jetz" , 
  #   option = "I", 
  #   direction = -1
  #       ) +
  scale_fill_carto_c(
    name = "jetz" ,
    type = "quantitative",
    palette = 18) + #18 19
  # scale_fill_gradient2(
  #   low = blue_grad[1], mid = blue_grad[2], high = blue_grad[3], 
  #   midpoint = 0.5
  # ) +
  
  coord_equal() +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) 



map_mpd <- 
sites_db %>% 
  ggplot() + 
  geom_raster(aes(x = x, y = y, fill = ses_mpd)) + 
  scale_fill_carto_c(
    name = "ses_mpd" ,
    type = "quantitative",
    palette = 5) + #18 19

  # scale_fill_gradient2(
  #   low = pink_grad[1], mid = pink_grad[2], high = pink_grad[3], 
  #   midpoint = -1.5
  # ) +
  coord_equal() +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

data_bi <- 
  sites_db %>% 
  bi_class(x = ses_mpd, y = prop_insitu, style = "equal", dim = 3)

map <- 
  ggplot() +
  geom_raster(data = data_bi, mapping = aes(x = x, y = y, fill = bi_class), show.legend = FALSE) +
  bi_scale_fill(pal = "GrPink2", dim = 3) +
  theme_bw() +
  coord_equal()

legend <- 
  bi_legend(pal = "GrPink2",
                    dim = 3,
                    ylab = "In Situ Div",
                    xlab = "SES_MPD ",
                    size = 8)

bi_map <- 
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.1, .1, 0.3, 0.3)


design <- 
"AB
 CC
 CC"

map_insitu_mpd <- 
wrap_plots(
  map_insitu, 
  map_mpd, 
  bi_map, 
  design = design
)

ggsave(
  "insitu.png", 
  map_insitu_mpd, 
  width = 7, 
  height = 11
)


sites_db %>% 
  ggplot(aes(x = mean_age_arrival,  y = prop_insitu)) +
  geom_point()

cor(sites_db$mean_age_arrival, sites_db$prop_insitu)
