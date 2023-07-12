
# Load packages -----------------------------------------------------------
{
library(tidyverse)
library(here)
library(terra)
library(sf)
library(ggridges)
library(ggdist)
library(gghalves)
library(patchwork)
library(glue)
}

# load data ---------------------------------------------------------------

env_files <- list.files(here("data", "env"), full.names = T)

altitute <- rast(env_files[1])
env_layers <- rast(env_files[-1])
env_layers_05 <- terra::aggregate(env_layers, 3, na.rm=TRUE)


evoregion_df <- read.csv(
  here("output", "evoregion", "evoregions_stand_names_df.csv")
  )


# prepare data ------------------------------------------------------------

alt_df <- terra::extract(altitute, evoregion_df[,1:2])
env_df <- terra::extract(env_layers_05, evoregion_df[,1:2]) 

names(alt_df) <- c("ID", "alt")
names(env_df) <- c("ID", "temp_mean", "prec_d_month", "prec_d_quarter", "temp_season")

evo_env_df <- 
bind_cols(evoregion_df, alt_df, env_df[,-1]) %>% 
  select(-ID)


# anova test for env difference among evoregions --------------------------

evo_test <- evo_env_df #%>% filter(phy_con %in% c("A", "B"))  

alt_test <- aov(log(evo_test$alt) ~ evo_test$phy_con) 
temp_test <- aov((evo_test$temp_mean/10) ~ evo_test$phy_con)
prec_test <- aov(log(evo_test$prec_d_quarter+1) ~ evo_test$phy_con)

# help function to add anova stats on plots
aov_stat <- function(aov_res)   {
  
  f_stat <- broom::tidy(aov_res)$statistic[1] |> 
    round(digits = 2) |>
    as.character()
  
  p_value <- broom::tidy(aov_res)$p.value[1] |>  
    round(digits = 4)
  
  if(p_value < 0.001) p_value <- "< 0.001"
  else p_value <- p_value |> as.character()
  
  c(f_stat = f_stat, p_value = p_value)
}


# |- tukey's test ----

# help function to add letters to plot 
# ref: https://r-graph-gallery.com/84-tukey-test.html
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(
    multcompView::multcompLetters(Tukey.levels, Letters = ".")['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# |-- altitude ----
tuk.alt <- TukeyHSD(alt_test, conf.level=.95)
tuk.alt.labels <- generate_label_df(tuk.alt , "evo_test$phy_con") |> 
  add_column(
    x = 1:5, y = 9
  )

# |-- temperature ----
tuk.temp <- TukeyHSD(temp_test, conf.level=.95)
tuk.temp.labels <- generate_label_df(tuk.temp , "evo_test$phy_con") |> 
  add_column(
    x = 1:5, y = 30
  )


# plots -------------------------------------------------------------------

theme_set(theme_bw())

colors_evo <-  c(
  "#363870",
  "#589eab",
  "#d3a838",
  "#362401",
  "#BB4455"
)

box_slab_plot <- list(
  geom_boxplot(
    size = 0.7, width = .2, alpha = .4,
    outlier.shape = 1
  ),
    # stat_slab(
    #   adjust = .7, ## bandwidth
    #   width = .7,
    #   alpha = 0.9,
    #   position = position_nudge(x = 0.2)
    # ),
  geom_density_ridges2(
    position = position_nudge(y = 0.2),
    scale = .7,
    alpha = 0.9,
    rel_min_height = 0.01
  ),
    scale_color_manual(values = colors_evo, guide = "none"),
    scale_fill_manual(values = colors_evo, guide = "none"),
  coord_flip(ylim = c(1,5.5))
)

# |-- altitude ----
(plot_alt <- 
ggplot(
  data = evo_env_df, # %>% filter(phy_con %in% c("A", "B")),
  aes(y = phy_con, x = log(alt), color = phy_con, fill = phy_con)
  ) +
  box_slab_plot +
  labs(
    title = "Altitude",
    subtitle = glue(
      "F-stat = {aov_stat(alt_test)[1]}\np-value {aov_stat(alt_test)[2]}"
      ),
    y = "Evoregions",
    x = "Log(Altitude [m])"
  ) +
    annotate(
      "text", 
      y = tuk.alt.labels$x, 
      x = tuk.alt.labels$y, 
      label = tuk.alt.labels$Letters, 
      size = 7
    ) +
  theme(
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 8)
  )
)


# |-- temperature ----
(plot_tmean <- 
  ggplot(
    evo_env_df, # %>% filter(phy_con %in% c("A", "B")),
    aes(y = phy_con, x = temp_mean/10, color = phy_con, fill = phy_con)
  ) +
  box_slab_plot +
  scale_x_continuous(n.breaks = 10) +
  labs(
    title = "Temperature",
    subtitle = glue(
      "F-stat = {aov_stat(temp_test)[1]}\np-value {aov_stat(temp_test)[2]}"
    ),
    y = "Evoregions",
    x = "Annual Mean Temperature [ºC]"
  )  +
    annotate(
      "text", 
      y = tuk.temp.labels$x, 
      x = tuk.temp.labels$y, 
      label = tuk.temp.labels$Letters, 
      size = 7
    ) +
  theme(
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 8)
  )
)




# save plot ---------------------------------------------------------------


env_chart <- 
plot_alt + plot_tmean +
  plot_annotation(tag_levels = "A")


ggsave(
  here("output", "fig", "04_evo_env_chart.png"),
  env_chart, 
  width = 9, 
  height = 6.5
)
  


















