
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(rcartocolor)

# source(here("function", "stand_evoregion_names.R"))
# source(here("function", "site_freq_evoreg.R"))

# load data ---------------------------------------------------------------

site_xy <- read.csv(here("data", "W_xy.csv"))

path_evo_consensus <- here(
  "output", "evoregion", "phy_consensus", "res_evoregion_phy_consensus.rds"
)

evo_consensus <- readRDS(path_evo_consensus)
evoregions_grp <- evo_consensus$Cluster_Evoregions

# path_evo_200 <- list.files(
#   here("output", "evoregion", "phy_200"), 
#   full.names = T
# )

# evoregions_df <- 
# map_dfc(c(path_evo_consensus, path_evo_200), function(x){
#   readRDS(x)$Cluster_Evoregions
# }) %>%  
#   set_names(
#     c("phy_con", paste0("phy_", 1:200))
#   ) %>% 
#   as.data.frame()


# standardize evoregion names ---------------------------------------------

# find the variability in number of evoregion among phylogenies
# n_evoregion <- map_int(evoregions_df, function(x){
#  levels(x)|> length()
# }) 
# 
# # there is a peak in five evoregion, and a maximum of 9 evoregion 
# hist(n_evoregion)

# |- defining the reference names for evoregions
pcps_vsc_thr <- evo_consensus$PCPS$vectors[,1:3]
k_grps_df <- data.frame(evo_grp = evoregions_grp, pcps_vsc_thr)

k_grps_df %>% 
  group_by(evo_grp) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  select(-evo_grp) %>% 
  dist() %>% 
  hclust("ward.D2") %>% 
  plot(hang = -1)

evo_names <- c("E", "C", "D", "B", "A")
evo_gr_names <- evo_names[evoregions_grp] %>%  as.factor()



blue_gold_red_2 <- c(
  "#363870",
  "#589eab",
  "#d3a838",
  "#362401",
  "#BB4455", #"#cb3b2e",
  "#f9c0c2"
)

evo_ref_df <- bind_cols(site_xy, phy_con = evo_gr_names)
evo_ref_df %>% 
  ggplot(aes(x, y, fill = phy_con)) +
  geom_raster(show.legend = T) +
  scale_fill_manual(
    values = blue_gold_red_2[c(2,1,3,5,4)]
  ) +
  labs(title = "Evoregions") +
  coord_equal() +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


# 
# evo_df_std <- stand_evoregion_names(evoregions_df, evo_ref)  %>%  
#   set_names(
#     c("phy_con", paste0("phy_", 1:200))
#   ) %>% 
#   as.data.frame()

evo_df_std_xy <- bind_cols(site_xy, phy_con = evo_gr_names)

write.csv(
  evo_df_std_xy, 
  here("output", "evoregion", "evoregions_stand_names_df.csv"), 
  row.names = FALSE
)

# evo_df_xy <- bind_cols(site_xy, evoregions_df)
# 
# not_five_evoregions <- 
# apply(evoregions_df, 2, function(x){
#   length(unique(x)) != 5
# })  %>% 
#   which() %>% 
#   names()
# 
# 
# 
# 
# #  |- check standardization -----------------------------------------------
# evo_df_std_xy <- read.csv(
#   here("output", "evoregion", "evoregions_stand_names_df.csv"), 
#   stringsAsFactors = T)
# # map theme
# 
# blue_gold_red_test <- c(
#   "#363870",
#   "#589eab",
#   "#d3a838",
#   "#362401",
#   "#cb3b2e",
#   "#f9c0c2",
#   "#346F3B",
#   "#8E8E8A",
#   "#1F1F1B"
# )
# 
# names(blue_gold_red_test) <- LETTERS[1:9]
# 
# greys <- c(
#   "#040400",
#   "#1F1F1B",
#   "#3B3B37",
#   "#575753",
#   "#73736F",
#   "#8E8E8A",
#   "#AAAAA6",
#   "#C6C6C2",
#   "#E2E2DE",
#   "#FEFEFA"
# )
# 
# my_theme <- list(
#   theme(
#     panel.background = element_rect(fill = "#FAF8F4"), 
#     panel.grid = element_blank(), 
#     text = element_text(color = greys[1]), 
#     title = element_text(color = greys[2]),
#     axis.text = element_text(color = greys[2]), 
#     axis.ticks = element_line(color = greys[3]), 
#     panel.border = element_rect(color = greys[4], fill = NA), 
#     axis.title = element_blank()
#   )
# )
# coast <- rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)
# 
# map.limits <- list(
#   x = c(-100, -30),
#   y = c(-35, 25)
# )
# 
# # map function
# evo_plot <- function(evo_data, col_name){
#   
#   ggplot() +
#     geom_tile(
#       data = evo_data, 
#       aes(x = x, y = y, fill = .data[[col_name]])
#     ) +
#     scale_fill_manual(
#       name = "Evoregion",
#       values = (blue_gold_red_test),
#       na.value = "#C6C6C2"
#     ) +
#     geom_sf(data = coast) +
#     coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
#     ggtitle(col_name) +
#     theme(
#       legend.position = "bottom"
#     ) +
#     my_theme
# }
# 
# ggplot() +
#   geom_tile(
#     data = evo_df_std_xy, 
#     aes(x = x, y = y, fill = phy_con)
#   )
# map_evo_ref <- evo_plot(evo_df_std_xy, "phy_con")
# 
# # |-- generate map comparisons ----
# 
# for(i in seq_along(not_five_evoregions)){
#   map_evo_compare <- evo_plot(evo_df_std_xy, not_five_evoregions[[i]])
#   Sys.sleep(2)
#   print(
#     map_evo_ref + map_evo_compare
#     )
# 
# }
# 
# 
# # frequence of classification to the same evoregion -----------------------
# 
# evo_df_std_xy <- read.csv(
#   here("output", "evoregion", "evoregions_stand_names_df.csv")
# )
# 
# freq_evo_reg <- site_freq_evoreg(evo_df_std_xy[,-c(1:3)])
# 
# evo_summary_df_xy <- bind_cols(site_xy, 
#                                uncertainty.idx = freq_evo_reg$uncertainty.idx)
# 
# map_evo_ref <- evo_plot(evo_df_std_xy, "phy_con")
# map_evo_freq <- evo_plot(evo_summary_df_xy, "uncertainty.idx")
# 
# map_evo_ref + map_evo_freq
# 
# 
# sum(!evo_df_std_xy$phy_con == evo_summary_df_xy$evo.var.stand)
# 
# #var_map <- 
# ggplot() +
#   geom_raster(
#     data = evo_summary_df_xy, 
#     aes(x = x, y = y, fill = .data[["uncertainty.idx"]])
#   ) +
#   scale_fill_carto_c(
#     name = "Uncertainty Index",
#     palette = "Temps",
#     limits = c(0, 1)
#   ) +
#   #scale_fill_manual(
#   #  name = "Evoregion",
#   #  values = (blue_gold_red_2[]),
#   #  na.value = "#C6C6C2"
#   #) +
#   geom_sf(data = coast) +
#   coord_sf(xlim = map.limits$x, ylim = map.limits$y) +
#   ggtitle("Phylogenetic Uncertainty in Evoregion Classification") +
#   theme(
#     legend.position = "bottom",
#     plot.title.position = "plot"
#   ) +
#   my_theme
# 
# 
#  var_map
# 
# 
