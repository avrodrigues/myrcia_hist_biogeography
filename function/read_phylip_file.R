
read_phylip_file <- function(path_to_phylip){
  phylip.mtx <- read.table(path_to_phylip, skip = 1, colClasses = rep("character", 2))
  phylip.area.names <- read.table(path_to_phylip, nrows = 1)
  
  rows <- phylip.area.names[,1]
  cols <- phylip.area.names[,2]
  
  areas <- phylip.area.names[,3:(2+cols)]  %>%  unlist()
  areas <- mapply(sub, pattern = "\\(|\\)", replacement = "", x = areas)
  names(areas) = NULL
  
  areas_01_ <- as.numeric(unlist(strsplit(phylip.mtx[,2], "")))  
  
  matrix(
    areas_01_, 
    rows,
    cols,
    byrow = T, 
    dimnames = list(
      phylip.mtx[, 1],
      areas
    )
  )
}


tipranges_from_phylip_tree <- function(path_to_phylip, tree){
  phylip_data <- read_phylip_file(path_to_phylip)
  
  areas <- colnames(phylip_data)
  tip_names <- rownames(phylip_data)
  
  tip_areas <- 
    sapply(1:nrow(phylip_data), function(i){
      row.areas <- areas[phylip_data[i,] == 1]
      paste0(row.areas, collapse = "")
    })
  
  tip_areas_order <- match(tree$tip.label, tip_names)
  tip_areas[tip_areas_order]
  
}


