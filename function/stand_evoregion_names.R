#' standardize evoregions names made from posterior sampling phylogenies
#' 
#' The function uses the reference classification to rename the classification 
#' of each phylogeny.  
#' 
#' 
#' @param evo.grp dataframe with evoregion original name for groups. In rows sites 
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
  
  if(max.evo != length(evo.names.ref)){
    stop(paste("Please provide", max.evo, "names in 'evo.names.ref'" ))
  }
  
  if(!inherits(evo.grp, "data.frame")){
    stop("'evp.grp' must be of class data.frame")
    }
  
  if(any(!sapply(evo.grp, inherits, "factor"))){
    stop("All columns in 'evo.grp' must be of class factor" )
  }
  
  phy_con <- names(evo.grp)[1]
  n_evo_con <- length(unique(evo.grp[[phy_con]]))
  
  
  evo_df_col_names <- names(evo.grp)
  names(evo_df_col_names) <- names(evo.grp)
  
  df_evo <- 
    map_dfc(evo_df_col_names, function(x){
      
      phy_otr <- x
      n_evo_otr <- length(unique(evo.grp[[phy_otr]]))
      
      # count the number of sites by combiantion of grps in 
      # phy_con and phy_otr
      evo_grp_1 <- 
        evo.grp %>% 
        select(.data[[phy_con]], .data[[phy_otr]]) %>% 
        count(.data[[phy_otr]], .data[[phy_con]]) %>% 
        arrange(.data[[phy_otr]], desc(n)) %>% 
        group_by(.data[[phy_otr]]) %>% 
        slice(1) %>%  
        # sort phy_con names by the number of sites in each combination 
        # (large n to small n)
        arrange(.data[[phy_con]], desc(n)) %>% 
        mutate(across(where(is.factor), as.integer))
      
      
      # if the number of grps in phy_otr are lower than in phy_con, 
      # then find duplicated rows and remove it
      if(n_evo_otr < n_evo_con){
        dup <- duplicated(evo_grp_1[[phy_otr]])
        
        
        # duplicated rows in phy_otr means that the group has sites in more than 
        # one group in phy_con. 
        # removing duplicated rows keep the grp in phy_otr with the same name as
        # the the group with more sites in phy_con
        evo_grp_1 <- evo_grp_1 %>% 
          filter(!duplicated(.data[[phy_otr]]))
      }
      else{
        dup <- duplicated(evo_grp_1[[phy_con]])
        
        
        if(any(dup)){
          n_ref_gr <- sum(!dup)
          w_dup <- which(dup)
          for(i in seq_along(w_dup))
            evo_grp_1[[phy_con]][w_dup[i]] <- n_ref_gr + i
        } 
      }
      
      new_names <- 
        evo_grp_1 %>% 
        mutate(new_names = evo.names.ref[.data[[phy_con]]]) %>% 
        arrange(.data[[phy_otr]]) %>% 
        pull(new_names)
      
      as.factor(new_names[evo.grp[[phy_otr]]])
    }) %>% 
    as.data.frame()
  
  return(df_evo)

}
  

 
  
 