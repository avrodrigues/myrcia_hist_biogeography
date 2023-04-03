
# load packages -----------------------------------------------------------

library(ape)
library(tidyverse)
library(here)


# load data ---------------------------------------------------------------

# species names
myrcia_binary_df <- read.csv(
  here("data", "distribution", "myrcia_binary_df_05_degree.csv")
)

names_distribution <- 
  myrcia_binary_df %>% 
  select(species) %>% 
  distinct() %>% 
  pull() 

# reviwed taxonomic names for phylogeny tips
pow.reviewed.df <- read.csv(
  here("data", "phylogeny", "phy_names_reviewed_pow_df.csv")
)

# trees
myrcia_tree_consensus <- read.tree(
  here("data", "phylogeny", "myrcia_group", "myrcia_clean.txt")
)

myrcia_tree_200 <- read.tree(
  here("data", "phylogeny", "myrcia_group", "myrcia_clean.trees.txt")
)


# Cleaning and taxonomic standadization -----------------------------------

##|- tips to remove ----

## remove all tips with taxa not identified at species level
tips_sp <- 
  pow.reviewed.df %>% 
  filter(is.na(name_accepted)) %>% 
  pull(phy_names)


## remove tips with duplicated revised accepted names
## - keep tip with the same name as the accepted vesion
## - keep tip without 'aff' to the accepted version
## - choose one tip randomly

accepted_duplicated <- 
pow.reviewed.df %>% 
  filter(
    duplicated(phy_new_names) | duplicated(phy_new_names, fromLast = TRUE)
  ) %>% 
  arrange(phy_new_names) %>% 
  pull(phy_new_names) %>% 
  unique()

remove.tip <- tibble(
  str_order = integer(),
  phy_names = character(),
  phy_new_names = character()
)

for(i in seq_along(accepted_duplicated)){
  temp_df <- 
  pow.reviewed.df %>% 
    filter(phy_new_names == accepted_duplicated[i]) %>% 
    select(str_order, phy_names, phy_new_names)
  
  ## - keep tip with the same name as the accepted vesion
  tip_accepted <- temp_df$phy_names == accepted_duplicated[i]
  if(any(tip_accepted)){
    remove.tip <- remove.tip %>% 
      add_row(temp_df[!tip_accepted, ])
    next()
  }
  
  tip_aff <- str_detect(temp_df$phy_names, "_aff_")
  if(any(tip_aff)){
    remove.tip <- remove.tip %>% 
      add_row(temp_df[tip_aff, ])
    next()
  }
  else{
    remove.tip <- remove.tip %>% 
      add_row(temp_df[1, ])
  }
}

# |-- object tips to remove ----
tips_to_remove <-
remove.tip %>% 
  pull(phy_names) %>% 
  c(tips_sp)

# |-- tips new names ----

tips_new_names_df <- 
  pow.reviewed.df %>% 
  filter(!phy_names %in% tips_to_remove)

# |- clean phylogenies
myrcia_tree_consensus$tip.label <- str_to_lower(myrcia_tree_consensus$tip.label )
myrcia_tree_consensus_cleaned <- drop.tip(
  myrcia_tree_consensus, 
  tips_to_remove
  )

myrcia_tree_200_cleaned <- vector("list", 200)

for(i in 1:200){
  
  myrcia_tree_200[[i]]$tip.label <- str_to_lower(
    myrcia_tree_200[[i]]$tip.label
    )
  
  myrcia_tree_200_cleaned[[i]] <- 
    drop.tip(
      myrcia_tree_200[[i]], 
      tips_to_remove
    )
}

# |- update taxonomic nomenclature

if(
  all(tips_new_names_df$phy_names == myrcia_tree_consensus_cleaned$tip.label)
){
  myrcia_tree_consensus_cleaned$tip.label <- tips_new_names_df$phy_new_names
}


for(i in 1:200){
  
  tips_df_order <- match(
    myrcia_tree_200_cleaned[[i]]$tip.label,
    tips_new_names_df$phy_names
  )
  
  if(
    all(tips_new_names_df$phy_names[tips_df_order] == 
        myrcia_tree_200_cleaned[[i]]$tip.label)
  ){
    myrcia_tree_200_cleaned[[i]]$tip.label <- 
      tips_new_names_df$phy_new_names[tips_df_order]
  }else{stop(paste("Error in", i))}
  
 
}


# macth tip names and species w/ distribution  ----------------------------


tips_without_dist <- 
tips_new_names_df$phy_new_names[
  !tips_new_names_df$phy_new_names %in% names_distribution
]



# save cleanned phylogenies -----------------------------------------------

dir_phy_cleaned <- here("data", "phylogeny", "phy_cleaned")

if(!dir.exists(dir_phy_cleaned)) dir.create(dir_phy_cleaned)


zero200 <- as.character(0:200 )
  
zero200_vec <- ifelse(
  nchar(zero200) == 1, 
  paste0("00", zero200), 
  ifelse(
    nchar(zero200) == 2, 
    paste0("0", zero200),
    zero200)
)

for(i in seq_along(zero200_vec)){
  if(i == 1){
    phy <- myrcia_tree_consensus_cleaned
    filename <- paste0(zero200_vec[i], "_phy_myrcia_cleaned_consensus.new")
  } else {
    phy <- myrcia_tree_200_cleaned[[i-1]]
    filename <- paste0(
      zero200_vec[i], "_phy_myrcia_cleaned_posterior_sampling.new"
      )
  }
  
  dir.save.phy <- here(dir_phy_cleaned, filename)
  phy.dropped <- drop.tip(phy, tips_without_dist)
  
  write.tree(phy.dropped, dir.save.phy)
}
