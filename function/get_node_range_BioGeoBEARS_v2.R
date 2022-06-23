#' Get node ranges from BioGeoBEARS biome reconstruction model
#'
#' @param BioGeoBEARS.result An object containing the result of BioGeoBEARS 
#'   reconstruction
#'   
#' @return Data frame with two columns, one indicating the node and other the ancestral biome range
#' 
#' @export

get_node_range_BioGeoBEARS_v2 <- function(BioGeoBEARS.result){
  
  
  geogfn <- BioGeoBEARS.result$inputs$geogfn
  max_range_size <- BioGeoBEARS.result$inputs$max_range_size
  include_null_range <- BioGeoBEARS.result$inputs$include_null_range
  
  relprobs_matrix <- BioGeoBEARS.result$ML_marginal_prob_each_state_at_branch_top_AT_node
  
  
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  areas = getareas_from_tipranges_object(tipranges)
  
  statenames = areas_list_to_states_list_new(areas, 
                                             maxareas = max_range_size, 
                                             include_null_range = include_null_range, 
                                             split_ABC = FALSE)
  
  MLstates = get_ML_states_from_relprobs(relprobs_matrix, 
                                         statenames, 
                                         returnwhat = "states")
  
  n_tips <- nrow(tipranges@df)
  
  anc_area <- MLstates[-c(1:n_tips)]
  node_names <- paste0("N", (n_tips+1):length(MLstates))
  
  data.frame(node_names, anc_area)
  
  
}
