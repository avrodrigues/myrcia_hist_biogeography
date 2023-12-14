# Load packages -----------------------------------------------------------
{
  library(tidyverse)
  library(here)
  library(terra)
  library(sf)
  library(spdep)
  library(nlme)
  library(glue)
  library(furrr)
}

# load data ---------------------------------------------------------------

env_files <- list.files(here("data", "env"), full.names = T)

altitute <- rast(env_files[1])
env_layers <- rast(env_files[-1])
env_layers_05 <- terra::aggregate(env_layers, 3, na.rm=TRUE)


evoregion_df <- read.csv(
  here("output", "evoregion", "evoregions_stand_names_df.csv")
)


# prepare data ------------------------------------------------------------

alt_df <- terra::extract(altitute, evoregion_df[,1:2])
env_df <- terra::extract(env_layers_05, evoregion_df[,1:2]) 

names(alt_df) <- c("ID", "alt")
names(env_df) <- c("ID", "temp_mean", "prec_d_month", "prec_d_quarter", "temp_season")

evo_test <- 
  bind_cols(evoregion_df, alt_df, env_df[,-1]) %>% 
  select(-ID)

evo_test_sf <- st_as_sf(evo_test, coords = c("x", "y"), crs = st_crs(4326))

ggplot(evo_test_sf) +
  geom_sf()

# create a weight matrix for distances -----

max_1nn <- max(unlist(nbdists(knn2nb(knearneigh(evo_test_sf, k=1)), evo_test_sf)))
nb_evo <- dnearneigh(evo_test_sf, d1 = 25, d2 = max_1nn)
listw_evo <- nb2listwdist(nb_evo, evo_test_sf, style = "W")

# Morans' I tests - anova residuals ----
anova_alt_test <- aov(log(evo_test$alt) ~ evo_test$phy_con) 
anova_temp_test <- aov((evo_test$temp_mean) ~ evo_test$phy_con)
anova_prec_test <- aov(log(evo_test$prec_d_quarter+1) ~ evo_test$phy_con)


resid_anova_alt_test <- residuals(anova_alt_test) %>% scale()
resid_anova_temp_test <- residuals(anova_temp_test) %>% scale()
resid_anova_prec_test <- residuals(anova_prec_test) %>% scale()


moransI_temp_anova <- moran.mc(resid_anova_temp_test, listw_evo, nsim = 999) 
moransI_prec_anova <- moran.mc(resid_anova_prec_test, listw_evo, nsim = 999)
moransI_alt_anova <- moran.mc(resid_anova_alt_test, listw_evo, nsim = 999) 


res_moran_df <- map(list(moransI_temp_anova, moransI_prec_anova, moransI_alt_anova), function(moran){
  broom::tidy(moran) %>% 
    mutate(model_name = word(moran$data.name, 1, sep = " "), .before = 1)
}) %>% list_rbind()

#removing spatial auto correlation (gls), then fit an ANOVA model -----
set.seed(2338)
evo_test_subset <- evo_test %>%
  slice_sample(prop = 0.1)


## find the best spatial correlation structure ----
cor_structure <- c("corExp", "corGaus", "corLin", "corRatio", "corSpher")
env_var <- c("temp_mean", "log(prec_d_quarter+1)", "log(alt)")

mod_eval <- function(mod){
  eval(parse(text = mod))
}

mod_gls <- mod_eval(glue("nlme::gls(temp_mean ~ phy_con, data = evo_test_subset, correlation = nlme::{cor_structure[1]}(form = ~x + y))"))

aov(mod_gls, data = evo_test_subset) |> summary()
qf(.95, 4, 465)

1-pf(18.22, 4, 465)

### temperature ----
# corSpher is not converging and was removed from the test


tictoc::tic()
plan(multisession, workers = 35)
list_cor_struc_temp <- future_map(1:10, function(x){
  
  evo_test_subset <- evo_test %>% 
    slice_sample(prop = 0.1)
  
  cor_structure <- c("corExp", "corGaus", "corLin", "corRatio", "corSpher")
  
  #cor_struc_temp <- 
  future_map(1:5, function(i) {
    model_call <-  glue("nlme::gls(temp_mean ~ phy_con, data = evo_test_subset, correlation = nlme::{cor_structure}(form = ~x + y))")
    mod_result <- tryCatch({mod_eval(model_call[i])}, error = function(e) {NULL})
    
    anova <- tryCatch({aov(mod_result, data = evo_test_subset) }, error = function(e) {NULL})
    if(!is.null(anova)){
      aov_res <- anova %>% summary()
      anova_df <- data.frame(F_stat = aov_res[[1]][1,"F value"], p_value = aov_res[[1]][1,"Pr(>F)"] ) 
    }else{
      anova_df <- data.frame(F_stat = NA, p_value = NA ) 
    }
    
    
    tukey <- tryCatch({TukeyHSD(anova)}, error = function(e) {NULL})
    if(!is.null(tukey)){
      tukey_diff <- tukey$phy_con[,"diff"]
      tukey_df <- data.frame(run = x, tukey_diff = tukey_diff, group = names(tukey_diff))
    }else{
      tukey_df <- data.frame(run = x, tukey_diff = NA, group = NA)
    }
    
    tibble(
      cor_structure = cor_structure[i],
      #model = list(mod_result),
      AIC = ifelse(is.null(mod_result), NA, AIC(mod_result)), 
      anova = list(anova_df), 
      tukey = list(tukey_df)
    )
    
  }) %>% 
    list_rbind() %>% 
    mutate(delta_AIC = AIC - min(AIC, na.rm = T))  %>% 
    arrange(delta_AIC) %>% 
    slice(1)
    #mutate(rank = row_number()) 
  
}, .options = furrr_options(seed = TRUE)) %>% list_rbind()
plan(sequential)
tocc <- tictoc::toc()


saveRDS(list_cor_struc_temp, "output/spat_anova_temp.rds")

## precipitation ----
plan(multisession, workers = 35)
list_cor_struc_prec <- future_map(1:1000, function(x){
  
  evo_test_subset <- evo_test %>% 
    slice_sample(prop = 0.1)
  
  cor_structure <- c("corExp", "corGaus", "corLin", "corRatio", "corSpher")
  env_var <- c("temp_mean", "log(prec_d_quarter+1)", "log(alt)")
  
  cor_struc_temp <- 
  future_map(1:5, function(i) {
    model_call <-  glue("nlme::gls(log(prec_d_quarter+1) ~ phy_con, data = evo_test_subset, correlation = nlme::{cor_structure}(form = ~x + y))")
    mod_result <- tryCatch({mod_eval(model_call[i])}, error = function(e) {NULL})
    
    
    anova <- tryCatch({aov(mod_result, data = evo_test_subset) }, error = function(e) {NULL})
    tukey <- tryCatch({TukeyHSD(anova)}, error = function(e) {NULL})
    
    tibble(
      cor_structure = cor_structure[i],
      model = list(mod_result),
      AIC = ifelse(is.null(mod_result), NA, AIC(mod_result)), 
      anova = list(anova), 
      tukey = list(tukey)
    )
    
  }) %>% 
    list_rbind() %>% 
    mutate(delta_AIC = AIC - min(AIC, na.rm = T))  %>% 
    arrange(delta_AIC) %>% 
    slice(1)
  #mutate(rank = row_number()) 
  
}, .options = furrr_options(seed = TRUE)) %>% list_rbind()
plan(sequential)

saveRDS(list_cor_struc_prec, "output/spat_anova_prec.rds")

## Altitude ----
plan(multisession, workers = 35)

list_cor_struc_alt <- future_map(1:1000, function(x){
  
  evo_test_subset <- evo_test %>% 
    slice_sample(prop = 0.1)
  
  cor_structure <- c("corExp", "corGaus", "corLin", "corRatio", "corSpher")
  env_var <- c("temp_mean", "log(prec_d_quarter+1)", "log(alt)")
  
  #cor_struc_temp <- 
  future_map(1:5, function(i) {
    model_call <-  glue("nlme::gls(log(alt) ~ phy_con, data = evo_test_subset, correlation = nlme::{cor_structure}(form = ~x + y))")
    mod_result <- tryCatch({mod_eval(model_call[i])}, error = function(e) {NULL})
    
    anova <- tryCatch({aov(mod_result, data = evo_test_subset) }, error = function(e) {NULL})
    tukey <- tryCatch({TukeyHSD(anova)}, error = function(e) {NULL})
    
    tibble(
      cor_structure = cor_structure[i],
      #model = list(mod_result),
      AIC = ifelse(is.null(mod_result), NA, AIC(mod_result)), 
      anova = list(anova), 
      tukey = list(tukey)
    )
    
  }) %>% 
    list_rbind() %>% 
    mutate(delta_AIC = AIC - min(AIC, na.rm = T))  %>% 
    arrange(delta_AIC) %>% 
    slice(1)
  #mutate(rank = row_number()) 
  
}, .options = furrr_options(seed = TRUE)) %>% list_rbind()
plan(sequential)


saveRDS(list_cor_struc_alt, "output/spat_anova_alt.rds")


