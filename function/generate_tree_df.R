# generate tree data frame ------------------------------------------------


generate_tree_df <- function(tree = NULL, n = NULL){
  
  
  require(ape)
  require(ggtree)
  require(treeio)
  
  
  if(is.null(tree)){
    if(is.null(n)) stop("provide a number of tips")
    tree <- ape::rtree(n)
    tree <- treeio::as.ultrametric(tree)
  }
  
  
  phy_plot <- ggtree(tree)
  
  phy_plot$data %>% 
    mutate(
      x = x-max(x),
      node_selected = FALSE,
      branch_selected = FALSE
    ) 
}

