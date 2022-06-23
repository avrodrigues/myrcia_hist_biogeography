interactive_phy <- function(phy){
  require(ggtree)
  require(plotly)
  
  gg_phy <- ggtree(phy) +
    geom_point(aes(text = paste("Node:", node)))
  ggplotly(gg_phy, tooltip = "text") %>%   
  layout(dragmode = "select")
}
