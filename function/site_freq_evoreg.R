#' Measure the phylogentic uncertainty of evoregion classification 
#' 
#' Given a set of classification for each site (with standard group names), the 
#' function caluculates the frequency of classification of a site in each 
#' evoregion. Then, measures an uncertainty index (see details) of the evoregion 
#' classification
#' 
#' @param site.evoreg a matrix or data frame with standardized group names. In 
#'   rows are the sites in columns are the different classifications. The group 
#'   in each column must to represent the same group. The function 
#'   `stand.evoregion.names` is intented to standardize the names across multiple 
#'   classification. 
#'   
#' @details The uncertainty index is calculated as `1 - (sd.site/sd.max)`. 
#'   `sd.site` is the standard deviation value for the frenquecy of classification 
#'   in each region for a site; `sd.max` is the maximum standard deviation value
#'   possible in a vector of a frequencies with length `n`. That is, to have a 
#'   maximum standard deviation value, one region hav a frequency of 1, and the 
#'   others have a value of 0. These indicates that all classifications assigned 
#'   the site to the same region. Therefore, `(sd.site/sd.max)` is a measure o 
#'   certainty of classification, then `1 - (sd.site/sd.max)` turn it in an 
#'   uncertainty value. The uncertainty index varies between 0 and 1. 
#'   
#'  @return a list of length two. $mtx.freq.evoreg is the matrix of frequencies 
#'  and $uncertainty.idx a numeric vector with the uncertainty values. 
#'   
#'  
site_freq_evoreg <- function(site.evoreg){
  # frequencia por sitio
  l.freq.evoreg <- lapply(1:nrow(site.evoreg), function(site){
    table(t(site.evoreg[site,]))/ncol(site.evoreg)
  })
  
  max.evoreg <- 
    sapply(site.evoreg, function(x){
      length(unique(x))
      })
  
  max.n.evoreg <- max(max.evoreg)
  whichmax.evoreg <- which.max(max.evoreg)
  
  mtx.freq.evoreg <- matrix(
    0, 
    nrow = nrow(site.evoreg), 
    ncol = max.n.evoreg
    )
  
  colnames(mtx.freq.evoreg) <- unique(site.evoreg[,whichmax.evoreg])
  row.names(mtx.freq.evoreg) <- 1:nrow(site.evoreg)
  
  for(i in 1:nrow(mtx.freq.evoreg)){
    temp.freq<- l.freq.evoreg[[i]]
    temp.col <- colnames(mtx.freq.evoreg) %in% names(temp.freq)
    
    mtx.freq.evoreg[i,temp.col] <- temp.freq
  }
  
  # incerteza padronizada (0,1) por sítio
  sd.freq <- apply(mtx.freq.evoreg, 1, sd)
  
  max.variation <- c(1, rep(0, ncol(mtx.freq.evoreg)-1))
  sd.max <- sd(max.variation)
  
  sd.stand <- 1-(sd.freq/sd.max)
  
  list(mtx.freq.evoreg = mtx.freq.evoreg, uncertainty.idx = sd.stand)
}
