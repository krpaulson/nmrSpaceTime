
# Fit NMR models

library(INLA)
library(tidyverse)
library(magrittr)
library(spdep)
library(stringr)
library(openxlsx)

set.seed(1528)


# settings ----------------------------------------------------------------

settings <- yaml::read_yaml("run_settings.yaml")
list2env(settings, .GlobalEnv)

years <- c(year_start_estimation:year_end_estimation)

survey_meta <- readRDS(paste0("Data/", country, "/metadata.rds"))
hold_out_year_start <- max(as.numeric(survey_meta$survey_year)) - 3
hold_out_years <- c(hold_out_year_start:2020)
# hold_out_area <- 1


# setup -------------------------------------------------------------------

# load GADM
dsn <- paste0("Data/", country, "/shapeFiles_gadm")
layer <- list.files(dsn)
layer <- layer[grep("_1.shp", layer)]
layer <- gsub(".shp", "", layer)
poly_adm1 <- rgdal::readOGR(dsn = dsn, layer = layer)

# create adjacency matrices for admin1 and admin2
admin1_mat <- poly2nb(SpatialPolygons(poly_adm1@polygons))
nb2INLA(file = paste0("Data/", country, "/admin1_", tolower(country_code), ".graph"), admin1_mat)
admin1_mat <- nb2mat(admin1_mat, zero.policy = TRUE)

# make row and column names NAME_1
colnames(admin1_mat) <- rownames(admin1_mat) <- poly_adm1$NAME_1
admin1_names <- data.frame(GADM = poly_adm1@data$NAME_1,
                           Internal = rownames(admin1_mat))

# load data
mod.dat <- readRDS(paste0("Data/", country, "/nmr_data_prepped.rds"))
start.year <- min(years)
mod.dat <- mod.dat %>%
  mutate(year = as.numeric(as.character(years))) %>% 
  filter(year >= start.year & year >= survey_year - 10)
if (is.numeric(mod.dat$urban)) {
  mod.dat$urban_indicator <- as.numeric(mod.dat$urban == 1)
  mod.dat$rural_indicator <- as.numeric(mod.dat$urban == 2)
} else {
  mod.dat$urban_indicator <- as.numeric(mod.dat$urban == "urban")
  mod.dat$rural_indicator <- as.numeric(mod.dat$urban == "rural")
}

# create admin1 variable which is numeric admin 1 alphabetically
mod.dat$admin1 <- as.numeric(factor(mod.dat$admin1.name))

# create binomial dataframe with counts by cluster
binom_df <- mod.dat %>%
  mutate(year = as.numeric(as.character(years))) %>%
  group_by(cluster, year, survey_year) %>%
  summarise(direct = sum(Y*v005)/sum(v005),
            Y = sum(Y),
            N = n(),
            admin1 = admin1 %>% nth(1),
            admin1.name = admin1.name %>% nth(1),
            admin1.char = admin1.char %>% nth(1),
            urban_indicator = urban_indicator %>% nth(1),
            rural_indicator = rural_indicator %>% nth(1)) %>%
  ungroup()

# expand to all location-times
all <- binom_df %>% tidyr::expand(admin1, year = years, urban_indicator, survey_year)
binom_df <- binom_df %>%
  dplyr::right_join(all)  %>%
  mutate(ay = paste0(stringr::str_pad(admin1, 2, pad = "0"), "_", year),
         rural_indicator = ifelse(urban_indicator == 1, 0, 1),
         survey_year = as.character(survey_year),
         year_copy = year) %>%
  as.data.frame()

# create spacetime.unstruct variable which is numeric ay alphabetically
binom_df <- binom_df %>% arrange(admin1, year)
binom_df$spacetime.unstruct <- as.numeric(factor(binom_df$ay))

# load urban proportions
urb_prop <- readRDS(paste0("Results/", country, "/admin1_urban_weights.rds"))
ref.tab <- read.xlsx(paste0("Data/", country, "/reference_table.xlsx"))
ref.tab <- ref.tab %>%
  dplyr::select(admin1.name = GADM, region = Internal) %>%
  unique()
urb_prop <- urb_prop %>%
  left_join(ref.tab, by = "region") %>%
  arrange(admin1.name, years)


# Type IV interaction -----------------------------------------------------

# code from Taylor Okonek

fname <- paste0("Data/", country, "/admin1_", tolower(country_code), ".graph")
graphfile <- read.delim2(file = fname)
icar_prec <- matrix(0, nrow = nrow(graphfile), ncol = nrow(graphfile))

for (i in 1:nrow(graphfile)) {
  vals <- graphfile[i,]
  vals <- vals %>% str_split(" ") %>% unlist 
  vals <- as.numeric(vals[3:length(vals)])
  icar_prec[i,vals] <- -1
}
diag(icar_prec) <- (apply(icar_prec, 2, sum) * -1)

# check that rows and columns sum to 0
apply(icar_prec, 1, sum)
apply(icar_prec, 2, sum)

# now get rw1 precision matrix
S <- length(years)
RW_prec <- INLA:::inla.rw(n = S, 
                          order = 1, 
                          scale.model = FALSE, # set scale.model  = F because we'll scale the interaction instead
                          sparse = TRUE)

# Kronecker product between ICAR x RW1
# order matters here! put the icar first if you want your spacetime random effects to be in
# order arrange(region, time)
R <- icar_prec %x% RW_prec

# get constraint matrix
eigen_things <- eigen(R)
# check that the appropriate number of eigen values are 0
length(which(eigen_things$values < 0.000001)) == (nrow(icar_prec) + S - 1)
tmp <- t(eigen_things$vectors[,(eigen_things$values < 0.000001)])
constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))

# scale R
R <- inla.scale.model(R, constr = constr.st)


# INLA model --------------------------------------------------------------

formulas <- list()
results <- list()
time_model_list <- c("rw1", "rw2", "ar1",
                     "rw1_iv", "rw2_iv", "ar1_iv")

formulas[["rw1"]] <- 
  Y ~ urban_indicator + rural_indicator - 1 +                   # urban/rural intercept
  f(survey_year, model = "iid", constr = TRUE) +                # random effect for survey year
  f(admin1, graph = fname, model = "bym2") +                    # bym2 space effect
  f(year, model = "rw1", constr = TRUE) +                       # rw1 time effect
  f(year_copy, model = "iid") +
  f(spacetime.unstruct, model = "iid")

formulas[["rw2"]] <- 
  Y ~ urban_indicator + rural_indicator - 1 +                   # urban/rural intercept
  f(survey_year, model = "iid", constr = TRUE) +                # random effect for survey year
  f(admin1, graph = fname, model = "bym2") +                    # bym2 space effect
  f(year, model = "rw2", constr = TRUE) +                       # rw2 time effect
  f(year_copy, model = "iid") +
  f(spacetime.unstruct, model = "iid")

formulas[["ar1"]] <- 
  Y ~ urban_indicator + rural_indicator - 1 +                   # urban/rural intercept
  f(survey_year, model = "iid", constr = TRUE) +                # random effect for survey year
  f(admin1, graph = fname, model = "bym2") +                    # bym2 space effect
  f(year, model = "ar1") +                                      # ar1 time effect
  f(spacetime.unstruct, model = "iid")

formulas[["rw1_iv"]] <- 
  Y ~ urban_indicator + rural_indicator - 1 +                   # urban/rural intercept
  f(survey_year, model = "iid", constr = TRUE) +                # random effect for survey year
  f(admin1, graph = fname, model = "bym2") +                    # bym2 space effect
  f(year, model = "rw1", constr = TRUE) +                       # rw2 time effect
  f(year_copy, model = "iid") +
  f(spacetime.unstruct,
    model = "generic0",
    Cmatrix = R,
    extraconstr = constr.st,
    rankdef = nrow(icar_prec) + S - 1)

formulas[["rw2_iv"]] <- 
  Y ~ urban_indicator + rural_indicator - 1 +                   # urban/rural intercept
  f(survey_year, model = "iid", constr = TRUE) +                # random effect for survey year
  f(admin1, graph = fname, model = "bym2") +                    # bym2 space effect
  f(year, model = "rw2", constr = TRUE) +                       # rw2 time effect
  f(year_copy, model = "iid") +
  f(spacetime.unstruct,
    model = "generic0",
    Cmatrix = R,
    extraconstr = constr.st,
    rankdef = nrow(icar_prec) + S - 1)

formulas[["ar1_iv"]] <- 
  Y ~ urban_indicator + rural_indicator - 1 +                   # urban/rural intercept
  f(survey_year, model = "iid", constr = TRUE) +                # random effect for survey year
  f(admin1, graph = fname, model = "bym2") +                    # bym2 space effect
  f(year, model = "ar1") +                                      # rw2 time effect
  f(spacetime.unstruct,
    model = "generic0",
    Cmatrix = R,
    extraconstr = constr.st,
    rankdef = nrow(icar_prec) + S - 1)

# if only one survey, remove survey effect
if (length(unique(binom_df$survey_year)) > 1) {
  survey_effect <- T
} else {
  survey_effect <- F
  for (tm in time_model_list) {
    tt <- terms(formulas[[tm]])
    survey_index <- which(attr(tt, "term.labels") %>% str_detect("survey"))
    tt <- drop.terms(tt, survey_index)
    formulas[[tm]] <- reformulate(attr(tt, "term.labels"), response = "Y", intercept = F)
  }
}

for (hold_out_area in 1:length(unique(binom_df$admin1))) {
  print(hold_out_area)

  # apply holdouts for validation
  binom_df_holdouts <- binom_df
  binom_df_holdouts[binom_df_holdouts$year %in% hold_out_years |
                      binom_df_holdouts$admin1 %in% hold_out_area, ]$Y <- NA
  
  # fit all models
  for (time_model in time_model_list) {
    print(time_model)
    results[[time_model]] <- inla(
      formula = formulas[[time_model]],
      Ntrials = N,
      data = binom_df_holdouts,
      family = "betabinomial",
      control.predictor = list(compute = TRUE),
      control.compute = list(return.marginals = TRUE,
                             config = TRUE,
                             return.marginals.predictor = TRUE))
  }
  
  
  # Samples -----------------------------------------------------------------
  
  expit_medians <- list()
  
  for (time_model in time_model_list) {
    
    print(time_model)
    
    result <- results[[time_model]]
    
    # obtain posterior samples
    nsamp <- 1000
    samp <- inla.posterior.sample(n = nsamp, result = result)
    
    # set up matrices to contain posterior samples for each term in our linear predictor
    # which we will later need to combine
    # only take the first half of the region id's since the bym2 returns c(total, spatial)
    region_idx <- which(rownames(samp[[1]]$latent) %>%
                          str_detect("admin1"))[1:length(unique(binom_df$admin1))]
    yearx <- which(rownames(samp[[1]]$latent) %>% str_detect("^year:"))
    if (time_model %like% "rw") yearx_iid <- which(rownames(samp[[1]]$latent) %>% str_detect("^year_copy"))
    ayx <- which(rownames(samp[[1]]$latent) %>% str_detect("spacetime.unstruct")) # interaction
    
    # urban intercept
    urbanx <- which(rownames(samp[[1]]$latent) %>% str_detect("urban"))
    ruralx <- which(rownames(samp[[1]]$latent) %>% str_detect("rural"))
    
    region_mat <- matrix(0, nrow = length(region_idx), ncol = nsamp)
    year_mat <- matrix(0, nrow = length(yearx), ncol = nsamp)
    if (time_model %like% "rw") year_iid_mat <- matrix(0, nrow = length(yearx_iid), ncol = nsamp)
    ay_mat <- matrix(0, nrow = length(ayx), ncol = nsamp)
    urban_mat <- matrix(0, nrow = 1, ncol = nsamp)
    rural_mat <- matrix(0, nrow = 1, ncol = nsamp)
    
    # fill in matrix with posterior samples
    for (i in 1:nsamp) {
      region_mat[,i] <- samp[[i]]$latent[region_idx]
      year_mat[,i] <- samp[[i]]$latent[yearx]
      if (time_model %like% "rw") year_iid_mat[,i] <- samp[[i]]$latent[yearx_iid]
      ay_mat[,i] <- samp[[i]]$latent[ayx]
      urban_mat[,i] <- samp[[i]]$latent[urbanx]
      rural_mat[,i] <- samp[[i]]$latent[ruralx]
    }
    
    # create space-time grid
    st_df <- expand.grid(region = 1:nrow(region_mat),
                         time = 1:nrow(year_mat)) %>%
      arrange(region, time)
    
    # prediction
    samp_mat_urban <-
      urban_mat[rep(1, nrow(ay_mat)),] +  # urban intercept
      region_mat[st_df$region,] +         # space effect
      year_mat[st_df$time,] +             # time effect
      ay_mat                              # space-time interaction
    samp_mat_rural <-
      rural_mat[rep(1, nrow(ay_mat)),] +  # rural intercept
      region_mat[st_df$region,] +         # space effect
      year_mat[st_df$time,] +             # time effect
      ay_mat                              # space-time interaction
    
    # iid time effect
    if (time_model %like% "rw") {
      samp_mat_urban <- samp_mat_urban + year_iid_mat[st_df$time,]
      samp_mat_rural <- samp_mat_rural + year_iid_mat[st_df$time,]
    }
    
    # now expit all of our samples to get back to p
    expit_samp_mat_urban <- SUMMER::expit(samp_mat_urban)
    expit_samp_mat_rural <- SUMMER::expit(samp_mat_rural)
    
    # use urban proportions to aggregate
    expit_samp_mat <-
      urb_prop$urban * expit_samp_mat_urban +
      urb_prop$rural * expit_samp_mat_rural
    
    # convert back to logit for SD on prediction scale
    # can use delta method later
    logit_samp_mat <- SUMMER::logit(expit_samp_mat)
    
    # compile summaries
    expit_medians[[time_model]] <- data.frame(
      time_model = time_model,
      nmr = apply(expit_samp_mat, 1, median),
      nmr_lower = apply(expit_samp_mat, 1, quantile, 0.025),
      nmr_upper = apply(expit_samp_mat, 1, quantile, 0.975),
      logit_nmr_sd = apply(logit_samp_mat, 1, sd),
      admin1_name = sort(unique(binom_df$admin1.name))[st_df$region],
      year = sort(unique(binom_df$year))[st_df$time]
    )
    
  }
  
  expit_medians <- do.call("rbind", expit_medians)
  
  # only save holdout
  hold_out_area_name <- sort(unique(binom_df$admin1.name))[hold_out_area]
  expit_medians <- expit_medians %>% filter(admin1_name == hold_out_area_name)
  
  saveRDS(expit_medians, paste0("Results/", country, "/results_holdout_", hold_out_area, ".rds"))

}
