select_clade <- function(tree_df, node_number){
  
  tree_df <- 
    tree_df %>% 
    mutate(
      node_selected = FALSE
    ) 
  
  if(is.null(node_number)){
    clade_df <- NULL
  }
  
  
  if(!is.null(node_number)){
    clade_df <- offspring(tree_df, node_number)
    
  }
  
  
  tree_df <- 
    tree_df %>% 
    mutate(
      node_selected = if_else(
        node %in% c(clade_df$node, node_number), 
        TRUE, 
        node_selected
      )
    )
  
  root_node <- sum(tree_df$isTip) + 1
  
  tree_df <- 
    tree_df %>% 
    mutate(
      branch_selected = map_chr(node, function(x){ 
        if(x == root_node) {
          return(tree_df %>% filter(node == root_node) %>% pull(node_selected))
        }
        parent(tree_df, x) %>% 
          pull(node_selected) %>%
          as.character()
      })
    )
  
  tree_df
}
