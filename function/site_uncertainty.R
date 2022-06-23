#' Calcualtes the uncetainty of evoregion classification for each site based in
#' multiple phylogenies
#'
#' @param evo_mtx a matrix with the evoregion classification for each site. Site 
#'   in rows, classification from each phylogeny in columns
#' 
#' @value a vector with the uncertainty

site_uncertainty <- function(evo_mtx){
  site_sim <- matrix(NA, nrow(evo_mtx), ncol(evo_mtx))
  site_uncertainty <- numeric(nrow(evo_mtx))
  
  
  for(i in 1:nrow(evo_mtx)){
    for(k in 1:ncol(evo_mtx)){
      site_evo <- evo_mtx[i,k]
      site_sim[,k] <- ifelse(site_evo == evo_mtx[,k], 1, 0)
    }
    site_sim
    site_uncertainty[i] <- apply(site_sim, 1, var) |> mean()
  }
  site_uncertainty
}