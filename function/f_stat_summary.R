f_stat_summary <- function(f_stat_df){
  require(glue)
  
  median_F <- f_stat_df$F_stat |> median()
  sd_F <- sd(f_stat_df$F_stat)
  p_value <- 1-pf(median_F, 4, 465)
  p_value <- ifelse(p_value < 0.001, "< 0.001", round(p_value,2))
  
  glue::glue("F-stat: {round(median_F,2)}±{round(sd_F,2)}, p-value: {(p_value)}")
}
