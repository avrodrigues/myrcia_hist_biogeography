#' standardize evoregions names made from posterior sampling phylogenies
#' 
#' The function uses the reference classification and to rename the classification 
#' of each phylogeny.  
#' 
#' 
#' @param evo.grp matrix with evoregion original name for groups. In rows sites 
#'   and in columns the different classification from each phylogeny. The first 
#'   column is used as reference for the standardization. 
#' @param evo.names.ref character vector with the name of evoregions for the 
#'   classification of reference, that is, the firt column of the \code{evo.grp}
#'   matrix.
#' 
#' 


stand_evoregion_names <- function(evo.grp, evo.names.ref){
  
  require(purrr)
  require(dplyr)
  
  max.evo <- max(sapply(evo.grp, function(x) length(unique(x))))
  
  if(max.evo > length(evo.names.ref)){
    stop(paste("Please provide at least", max.evo, "names in 'evo.names.ref'" ))
  }
  
  phy_con <- names(evo.grp)[1]
  
  df_evo <- 
    map_dfc(names(evo.grp), function(x){
      
      phy_otr <- x
      
      evo_grp_1 <- 
        evo.grp %>% 
        select(.data[[phy_con]], .data[[phy_otr]]) %>% 
        count(.data[[phy_otr]], .data[[phy_con]]) %>% 
        arrange(.data[[phy_otr]], desc(n)) %>% 
        group_by(.data[[phy_otr]]) %>% 
        slice_max(n) %>% 
        arrange(.data[[phy_con]], desc(n)) %>% 
        mutate(across(where(is.factor), as.integer))
      
      dup <- duplicated(evo_grp_1[[phy_con]])
      
      
      if(any(dup)){
        n_ref_gr <- sum(!dup)
        w_dup <- which(dup)
        for(i in seq_along(w_dup))
          evo_grp_1[[phy_con]][w_dup[i]] <- n_ref_gr + i
      } 
      
      new_names <- 
        evo_grp_1 %>% 
        mutate(new_names = evo_ref[phy_con]) %>% 
        arrange(.data[[phy_otr]]) %>% 
        pull(new_names)
      
      as.factor(new_names[evo.grp[[phy_otr]]])
    })
  
  return(df_evo)

}
  

 
  
 