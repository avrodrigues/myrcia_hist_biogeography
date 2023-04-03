library(Herodotools)
library(picante)


resDECJ <- readRDS(
  "output/biogeobears/list_results_models_biogeobears_phy_consenso.rds"
)$resDECj

resBayAreaJ <- readRDS(
  "output/biogeobears/list_results_models_biogeobears_phy_consenso.rds"
  )$resBAYAREALIKEj


# load data ---------------------------------------------------------------


phy.path <- here(
  "data", "phylogeny", "phy_cleaned", "000_phy_myrcia_cleaned_consensus.new"
)
myrcia_tree <- read.tree(phy.path)

geog.path <- here(
  "output", "biogeobears", "spp_area", "000_areas_myrcia_phy_consensus.data"
)

# evoregion classification
evo_df <- read.csv(here("output", "evoregion", "evoregions_stand_names_df.csv"))
site_xy <- evo_df[,1:2]
evo_mtx <- data.frame(evoregion = evo_df[, -c(1:2)])

# species composition
myrcia_comp <- read.csv(here("data", "W.csv"))


# assemblage age and in situ dibersification ---------------------------------

node.area <- 
  get_node_range_BioGeoBEARS(
    resDECJ,
    phyllip.file = geog.path,
    myrcia_tree,
    max.range.size = 2
  )

# calculating age arrival 
age_comm <- age_arrival(W = myrcia_comp, 
                        tree = myrcia_tree, 
                        ancestral.area = node.area, 
                        biogeo = evo_mtx) 

# calculating in situ diversification

myrcia_diversification <- 
  db_diversification(W = myrcia_comp,
                     tree = myrcia_tree, 
                     ancestral.area = node.area, 
                     biogeo = evo_mtx, 
                     diversification = "jetz",
                     type = "equal.splits")

# Calculating SES MPD
dis <- cophenetic(myrcia_tree) #cophenetic distance matrix
org<- SYNCSA::organize.syncsa(comm=myrcia_comp, phylodist=dis) #organizing matrices

ses_mpd_myrcia <- ses.mpd(
  org$community, dis = org$phylodist, null.model = "taxa.labels", nperm = 1000
)

mpd.obs.z <- ses_mpd_myrcia$mpd.obs.z

sites_db <- bind_cols(
  evo_df,
  age = age_comm$mean_age_per_assemblage,
  diversification_model_based = myrcia_diversification$model_based_Jetz_harmonic_mean_site,
  diversification = myrcia_diversification$Jetz_harmonic_mean_site
) %>% 
  mutate(
    prop_insitu = diversification_model_based/diversification,
    #ses_mpd = mpd.obs.z
  )

saveRDS(sites_db, "sites_db_metrics.rds")

