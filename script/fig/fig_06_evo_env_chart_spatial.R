
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

list_cor_struc_temp <- readRDS("output/spat_anova_temp.rds")
list_cor_struc_prec <- readRDS("output/spat_anova_prec.rds")
list_cor_struc_alt  <- readRDS("output/spat_anova_alt.rds")

# prepare data ------------------------------------------------------------

alt_df <- terra::extract(altitute, evoregion_df[,1:2])
env_df <- terra::extract(env_layers_05, evoregion_df[,1:2])

names(alt_df) <- c("ID", "alt")
names(env_df) <- c("ID", "temp_mean", "prec_d_month", "prec_d_quarter", "temp_season")

evo_env_df <-
bind_cols(evoregion_df, alt_df, env_df[,-1]) %>%
  select(-ID)
# 
# evo_test <- evo_env_df 
# 
# set.seed(2338)
# evo_test_subset <- evo_test |> 
#   slice_sample(prop = 0.1)

# anova test for env difference among evoregions --------------------------


# help function to add anova stats on plots
aov_stat <- function(f_stat_df)   {
  
  median_F <- f_stat_df$F_stat |> median()
  sd_F <- sd(f_stat_df$F_stat)
  p_value <- 1-pf(median_F, 4, 465)
  p_value <- ifelse(p_value < 0.001, "< 0.001", round(p_value,2))
  
  glue::glue("F-stat: {round(median_F,2)} ±{round(sd_F,2)}, p-value: {(p_value)}")
}






# |- tukey's test ----

tukey_test_df <- function(tukey_res_df){
  require(dplyr)
  
  tukey_res_df |> 
    group_by(group) |> 
    summarise(
      mean_diff = mean(tukey_diff),
      lower_quant = quantile(tukey_diff, probs = 0.025),
      upper_quant = quantile(tukey_diff, probs = 0.975),
      sig_diff = !(lower_quant < 0 & upper_quant > 0)
    )
}

# help function to add letters to plot 
# ref: https://r-graph-gallery.com/84-tukey-test.html
generate_label_df <- function(TUKEY_res_df){
  
  # Extract labels and factor levels from Tukey post-hoc 
  temp_res_t <- tukey_test_df(TUKEY_res_df)
  
  vec_diff <- temp_res_t$sig_diff
  names(vec_diff) <- temp_res_t$group
  
 
  Tukey.labels <- data.frame(
    multcompView::multcompLetters(vec_diff, Letters = letters)['Letters']
    )
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}



# |-- temperature ----
tuk.temp <- list_cor_struc_temp$tukey |>  list_rbind()
tuk.temp.labels <- generate_label_df(tuk.temp) |> 
  mutate(
    Letters = str_replace_all(
      Letters, c("a" = "●", "b" = "▲", "c" = "■")  
      ),
    x = 1:5,
    y = 30)

# |-- precipitation ----
tuk.prec <- list_cor_struc_prec$tukey |>  list_rbind()
tuk.prec.labels <- generate_label_df(tuk.prec) |> 
  mutate(
    Letters = str_replace_all(
      Letters, c("a" = "●", "b" = "▲", "c" = "■", "d" = "▼")  
    ),
    x = 1:5,
    y = 3000
  )

# |-- altitude ----
tuk.alt <- list_cor_struc_alt$tukey |>  list_rbind()
tuk.alt.labels <- generate_label_df(tuk.alt) |> 
  mutate(
    Letters = str_replace_all(
      Letters, c("a" = "●", "b" = "▲", "c" = "■")  
    ),
    x = 1:5,
    y = 8000
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
    scale = .6,
    alpha = 0.9,
    rel_min_height = 0.01
  ),
    scale_color_manual(values = colors_evo, guide = "none"),
    scale_fill_manual(values = colors_evo, guide = "none"),
  coord_flip(ylim = c(1,5.5)),
  theme(
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 12)
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
    title = "Annual Mean Temperature",
    subtitle = aov_stat(list_cor_struc_temp$anova |> list_rbind()),
    y = "Evoregions",
    x = "Temperature (ºC)"
  )  +
    annotate(
      "text", 
      y = tuk.temp.labels$x, 
      x = tuk.temp.labels$y, 
      label = tuk.temp.labels$Letters, 
      size = 5
    ) 
)

# |-- precipitation ----
(plot_prec <- 
   ggplot(
     evo_env_df, # %>% filter(phy_con %in% c("A", "B")),
     aes(y = phy_con, x = prec_d_quarter+1, color = phy_con, fill = phy_con)
   ) +
   box_slab_plot +
   labs(
     title = "Precipitaion of the Driest Quarter",
     subtitle = aov_stat(list_cor_struc_prec$anova |> list_rbind()),
     y = "Evoregions",
     x = "Precipitation + 1 (mm)"
   ) +
   scale_x_continuous(
     trans= "log",
     breaks = c(0,2,5,10,20,50,100,200,500,1000,2000)
   ) +
   annotate(
     "text", 
     y = tuk.prec.labels$x, 
     x = tuk.prec.labels$y, 
     label = tuk.prec.labels$Letters, 
     size = 5
   )  +
   theme(
     panel.grid.minor = element_blank()
   )
)


# |-- altitude ----
(plot_alt <- 
   ggplot(
     data = evo_env_df, # %>% filter(phy_con %in% c("A", "B")),
     aes(y = phy_con, x = alt, color = phy_con, fill = phy_con)
   ) +
   box_slab_plot +
   labs(
     title = "Altitude",
     subtitle = aov_stat(list_cor_struc_alt$anova |> list_rbind()),
     y = "Evoregions",
     x = "Altitude (Km)"
   ) +
   scale_x_continuous(
     trans= "log", 
     breaks = c(0, 10, 25, 100, 250, 500, 1000, 2500, 5000)
     ) +
   annotate(
     "text", 
     y = tuk.alt.labels$x, 
     x = tuk.alt.labels$y, 
     label = tuk.alt.labels$Letters, 
     size = 5
   ) +
   theme(
     panel.grid.minor = element_blank()
   )
) 



# save plot ---------------------------------------------------------------
env_chart <- 
plot_tmean + plot_prec + plot_alt +
  plot_annotation(tag_levels = "A")


ggsave(
  here("output", "fig", "05_evo_env_chart_spatial.png"),
  env_chart, 
  width = 12, 
  height = 5
)
  


















