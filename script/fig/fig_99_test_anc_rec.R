library(ggforce)
library(ggplot2)
library(dplyr)
library(cowplot)


evo_cols <-   c( 
  "A" = "#37396C",
  "B" = "#58A0A9", 
  "C" = "#D3A838",
  "D" = "#362401",
  "E" = "#BB4455")

area <- "ABC"

evo_diagram <- function(range){

  range_vec <- strsplit(range, "") %>% unlist()
  
  x <- 1:nchar(range)
  y <- rep(1, nchar(range))
  
  xlim <- c(mean(x)-2.5, mean(x)+2.5)
  
  df_tiles <- data.frame(x, y, evo = range_vec)
  
  diagram <- 
  ggplot(df_tiles, aes(x, y, fill = evo)) +
    geom_tile(show.legend = F) +
    geom_text(
      aes(x, y, label = evo), 
      color = "white"
    ) +
    scale_fill_manual(values = evo_cols) +
    coord_equal(xlim = xlim) +
    theme_void() +
    theme(
      plot.background = element_rect(colour = NA, fill = "transparent")
    )
  
  diagram
  
}

evo_diagram("BCE")

test_diagrams <- 
  map(tree_df$area, function(x){
  evo_diagram(x)
})

names(test_diagrams) <- 1:nrow(tree_df)

p <-
  ggtree(
    tree_df, 
    layout = "rect",
    aes(x = x), 
    size = 2,
    color = "#8E8E8A"
  )  +
  geom_tiplab(
    color = "#202015", 
    offset = 0.5, 
    fontface = 3
  ) +
  scale_x_continuous(
    limits = c(-30, 5), 
    breaks = seq(-30, 5, 5), 
    labels = c(paste0(seq(30, 0, -5), "Ma"), "") 
  ) +
  theme_tree2() 



new_chart_phy <- inset(p, test_diagrams, width=.09, height=.09)
ggsave(
  "output/fig/new_chart_phy.png",
  new_chart_phy , 
  width = 10, 
  height = 16
) 


