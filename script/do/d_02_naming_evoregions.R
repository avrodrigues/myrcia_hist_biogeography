
# load packages -----------------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(rcartocolor)

# load data ---------------------------------------------------------------

site_xy <- read.csv(here("data", "W_xy.csv"))

path_evo_consensus <- here(
  "output", "evoregion", "phy_consensus", "res_evoregion_phy_consensus.rds"
)

evo_consensus <- readRDS(path_evo_consensus)
evoregions_grp <- evo_consensus$cluster_evoregions

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

evo_names <- c("D", "B", "E", "A", "C")
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
    values = blue_gold_red_2[]
  ) +
  labs(title = "Evoregions") +
  coord_equal() +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


evo_df_std_xy <- bind_cols(site_xy, phy_con = evo_gr_names)

write.csv(
  evo_df_std_xy, 
  here("output", "evoregion", "evoregions_stand_names_df.csv"), 
  row.names = FALSE
)

