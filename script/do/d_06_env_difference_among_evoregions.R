# Load packages -----------------------------------------------------------
{
  library(tidyverse)
  library(here)
  library(terra)
  library(sf)
  library(spdep)
  library(nlme)
  library(glue)
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


resid_anova_alt_test <- residuals(anova_alt_test) |> scale()
resid_anova_temp_test <- residuals(anova_temp_test) |> scale()
resid_anova_prec_test <- residuals(anova_prec_test) |> scale()


moransI_temp_anova <- moran.mc(resid_anova_temp_test, listw_evo, nsim = 999)
moransI_prec_anova <- moran.mc(resid_anova_prec_test, listw_evo, nsim = 999)
moransI_alt_anova <- moran.mc(resid_anova_alt_test, listw_evo, nsim = 999) 

# removing spatial auto correlation (gls), then fit an ANOVA model -----
set.seed(2338)
evo_test_subset <- evo_test |> 
  slice_sample(prop = 0.1)


## find the best spatial correlation structure ----
cor_structure <- c("corExp", "corGaus", "corLin", "corRatio", "corSpher")
env_var <- c("temp_mean", "log(prec_d_quarter+1)", "log(alt)")

mod_eval <- function(mod){
  eval(parse(text = mod))
}

### temperature ----
# corSpher is not converging and was removed from the test
cor_struc_temp <- map(1:4, function(i) {
  model_call <-  glue("gls({env_var[1]} ~ 1, data = evo_test_subset, correlation = {cor_structure}(form = ~x + y))")
  mod_result <- mod_eval(model_call[i])
  
  data.frame(cor_structure = cor_structure[i], AIC = AIC(mod_result))
})  |> 
  list_rbind() |> 
  mutate(delta_AIC = AIC - min(AIC))  |> 
  arrange(delta_AIC)


### precipitation ----
cor_struc_prec <- map(c(1, 2, 4), function(i) {
  model_call <-  glue("gls({env_var[2]} ~ 1, data = evo_test_subset, correlation = {cor_structure}(form = ~x + y))")
  mod_result <- mod_eval(model_call[i])
  
  data.frame(cor_structure = cor_structure[i], AIC = AIC(mod_result))
})  |> 
  list_rbind() |> 
  mutate(delta_AIC = AIC - min(AIC))  |> 
  arrange(delta_AIC)

### altitude ----
cor_struc_alt <- map(1:5, function(i) {
  model_call <-  glue("gls({env_var[3]} ~ 1, data = evo_test_subset, correlation = {cor_structure}(form = ~x + y))")
  mod_result <- mod_eval(model_call[i])
  
  data.frame(cor_structure = cor_structure[i], AIC = AIC(mod_result))
})  |> 
  list_rbind() |> 
  mutate(delta_AIC = AIC - min(AIC))  |> 
  arrange(delta_AIC)

cor_struc_temp; cor_struc_prec; cor_struc_alt

## corExp is the optimal spatial correlation structure for all environmental variables ##

gls_temp_null <- gls(temp_mean ~ 1, data = evo_test_subset,
                     correlation = corExp(form = ~x + y))

gls_prec_null <- gls(log(prec_d_quarter+1) ~ 1, data = evo_test_subset,
                     correlation = corExp(form = ~x + y))

gls_alt_null <- gls(log(alt) ~ 1, data = evo_test_subset,
                    correlation = corExp(form = ~x + y))



resid_gls_temp_null <- gls_temp_null$residuals |> scale()
resid_gls_prec_null <- gls_prec_null$residuals |> scale()
resid_gls_alt_null <- gls_alt_null$residuals |> scale()


anova_temp_spat <- aov(resid_gls_temp_null ~ evo_test_subset$phy_con)
anova_prec_spat <- aov(resid_gls_prec_null ~ evo_test_subset$phy_con)
anova_alt_spat <- aov(resid_gls_alt_null ~ evo_test_subset$phy_con)

tuk_temp <- TukeyHSD(anova_temp_spat, conf.level=.95)
tuk_prec <- TukeyHSD(anova_prec_spat, conf.level=.95)
tuk_alt <- TukeyHSD(anova_alt_spat, conf.level=.95)

list_anova_results <- list(
  anova = list(temp = anova_temp_spat, prec = anova_prec_spat, alt = anova_alt_spat),
  tukey = list(temp = tuk_temp, prec = tuk_prec, alt = tuk_alt)
)

saveRDS(list_anova_results, "output/list_anova_results.rds")


