
# load packages -----------------------------------------------------------

library(tidyverse)
library(ggdist)
source("function/f_stat_summary.R")

# load data ---------------------------------------------------------------

list_cor_struc_temp <- readRDS("output/spat_anova_temp.rds")
list_cor_struc_prec <- readRDS("output/spat_anova_prec.rds")
list_cor_struc_alt  <- readRDS("output/spat_anova_alt.rds")


# prepare data ------------------------------------------------------------

list_cor_struc_temp$cor_structure |> table()

anova_temp_df <- list_cor_struc_temp$anova |> list_rbind()

f_stat_summary(anova_temp_df)

anova_temp_df |> 
  ggplot() +
  geom_histogram(aes(x = F_stat), bins = 20)

tukey_temp_df <- list_cor_struc_temp$tukey |> 
  list_rbind() |> 
  mutate(
    gr1 = word(group, 1, sep = "-"),
    gr2 = word(group, 2, sep = "-")
  )

ggplot(tukey_temp_df, aes(x = tukey_diff, y = group)) +
  stat_pointinterval(
    # aes(
    #   fill = after_stat(cut_cdf_qi(cdf, .width = c(0.95, 1)))),
    # show.legend = FALSE
    ) +
  #scale_fill_manual(values = c("lightblue", "gray50")) +
  geom_vline(xintercept = 0, linetype = 2) +
  #facet_wrap("group") +
  theme_bw()

# precipitation ----------------------
list_cor_struc_prec$cor_structure |> table()


anova_prec_df <- list_cor_struc_prec$anova |> list_rbind()
f_stat_summary(anova_prec_df)

anova_prec_df |> 
  ggplot() +
  geom_histogram(aes(x = F_stat), bins = 20)

tukey_prec_df <- list_cor_struc_prec$tukey |> 
  list_rbind() |> 
  mutate(
    gr1 = word(group, 1, sep = "-"),
    gr2 = word(group, 2, sep = "-")
  )

ggplot(tukey_prec_df, aes(x = tukey_diff, y = group)) +
  stat_pointinterval(
    # aes(
    #   fill = after_stat(cut_cdf_qi(cdf, .width = c(0.95, 1)))),
    # show.legend = FALSE
  ) +
  #scale_fill_manual(values = c("lightblue", "gray50")) +
  geom_vline(xintercept = 0, linetype = 2) +
  #facet_wrap("group") +
  theme_bw()

ggplot(tukey_prec_df, aes(x = tukey_diff)) +
  stat_halfeye(
    aes(
      fill = after_stat(cut_cdf_qi(cdf, .width = c(0.9, 1)))),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("lightblue", "gray70")) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_grid(gr1 ~ gr2) +
  theme_bw()

# altitude ----------------------
list_cor_struc_alt$cor_structure |> table()


anova_alt_df <- list_cor_struc_alt$anova |> list_rbind()
f_stat_summary(anova_alt_df)
anova_alt_df |> 
  ggplot((aes(x = F_stat)) +
  stat_pointinterval)

tukey_alt_df <- list_cor_struc_alt$tukey |> 
  list_rbind() |> 
  mutate(
    gr1 = word(group, 1, sep = "-"),
    gr2 = word(group, 2, sep = "-")
  )

ggplot(tukey_alt_df, aes(x = tukey_diff, y = group)) +
  stat_pointinterval(
    # aes(
    #   fill = after_stat(cut_cdf_qi(cdf, .width = c(0.95, 1)))),
    # show.legend = FALSE
  ) +
  #scale_fill_manual(values = c("lightblue", "gray50")) +
  geom_vline(xintercept = 0, linetype = 2) +
  #facet_wrap("group") +
  theme_bw()
