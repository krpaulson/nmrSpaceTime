
# Prep inputs

library(rdhs)
library(tidyverse)
library(SUMMER)
library(spdep)
library(geosphere)

source("helperFunctions.R")


# Settings ----------------------------------------------------------------

countries <- "Liberia"
country_codes <- "LBR"
year_start_country <- 2013 # first year's survey to include for this country, so that we only use 1 sampling frame
year_start <- 2000
year_end <- 2020

## Fill in your email and DHS project name if not yet configured
# set_rdhs_config(
#   email = EMAIL,
#   project = PROJECT_NAME
# )


# DHS data ----------------------------------------------------------------

# get country IDs
ids <- dhs_countries(returnFields=c("CountryName", "DHS_CountryCode"))
ids <- ids %>% filter(CountryName %in% countries)

# find all the surveys that fit our search criteria
survs <- dhs_surveys(countryIds = ids$DHS_CountryCode,
                     surveyYearStart = year_start)

# find datasets from these surveys
datasets <- dhs_datasets(surveyIds = survs$SurveyId)

# format
datasets <- datasets %>% 
  filter(FileType == "Geographic Data" | (FileType == "Births Recode" & FileFormat == "Stata dataset (.dta)")) %>%
  select(SurveyId, DatasetType, country = CountryName, survey_year = SurveyYear, file_name = FileName)

# format wide (one line per survey with GPS file and survey file)
datasets_wide <- datasets %>%
  mutate(file_name = sub("\\..*", "", file_name),
         DatasetType = sub(" Datasets", "", DatasetType)) %>%
  pivot_wider(names_from = DatasetType, values_from = file_name) %>%
  filter(!is.na(Survey))

# remove GPS if missing BR from datasets
datasets <- datasets %>% filter(SurveyId %in% datasets_wide$SurveyId)

# save metadata
saveRDS(datasets_wide, file = "Data/metadata.rds")

# download
downloads <- get_datasets(datasets$file_name, output_dir_root = "Data")

## NOTE: can start from here if data are already downloaded

# get survey metadata and subset to those using most-recent sampling frame
survey_meta <- readRDS("Data/metadata.rds")
survey_subset <- data.frame(country = countries, yr = year_start_country)
survey_meta <- survey_meta %>%
  right_join(survey_subset, by = "country") %>%
  filter(survey_year >= yr)

# loop through survey location-year and load data
mod_dat <- data.frame()
for (i in 1:nrow(survey_meta)) {
  
  print(paste(survey_meta[i,]$country, survey_meta[i,]$survey_year))
  
  # get survey data and format with SUMMER::get_births
  data.raw <- readRDS(paste0("Data/", survey_meta$Survey[i], ".rds"))
  
  # make sure alive variable is coded correctly for getBirths
  alive <- attr(data.raw$b5, which = "labels")
  names(alive) <- tolower(names(alive))
  data.raw$b5 <- ifelse(data.raw$b5 == alive["yes"][[1]], "yes", "no")
  data.raw$b5 <- factor(data.raw$b5, levels = c("yes", "no"))
  
  data <- suppressMessages(getBirths(
    data = data.raw,
    surveyyear = as.numeric(survey_meta$survey_year[i]),
    year.cut = seq(year_start, as.numeric(survey_meta$survey_year[i]) + 1, 1),
    strata = "v022",
    compact = F
  ))
  data$country <- survey_meta[i,]$country
  data <- data[ , c("country", "survey_year", "v001", "time", "age", "v005", "v025", "strata", "died")]
  data <- data[data$age == "0", ]
  
  # get GPS data if it exists, and use it to merge on location metadata
  gps_file <- paste0("Data/", survey_meta$GPS[i], ".rds")
  if (file.exists(gps_file)) {
    
    gps <- readRDS(gps_file)
    
    # get GADM shapefiles for Admin 1 and Admin 2
    gadm_dir <- paste0("Data/", survey_meta$country[i], "/shapeFiles_gadm")
    gadm_files <- list.files(gadm_dir)
    gadm1_file <- gsub(".shp", "", gadm_files[grep("_1.shp", gadm_files)])
    gadm2_file <- gsub(".shp", "", gadm_files[grep("_2.shp", gadm_files)])
    gadm1 <- rgdal::readOGR(dsn = gadm_dir, layer = gadm1_file, encoding = "UTF-8", use_iconv = TRUE, verbose = F)
    gadm2 <- rgdal::readOGR(dsn = gadm_dir, layer = gadm2_file, encoding = "UTF-8", use_iconv = TRUE, verbose = F)
    
    # use custom function to merge on all location metadata
    data <- add_adm1_adm2_meta(data, gps, gadm1, gadm2)
    
  } else {
    print("Missing GPS dataset")
    data[, c("LATNUM", "LONGNUM", "coords.x1", "coords.x2", "admin2", "admin2.char",
             "admin2.name", "admin1", "admin1.char", "admin1.name")] <- NA
  }
  
  data <- labelled::remove_labels(data, keep_var_label = T)
  mod_dat <- rbind(mod_dat, data)
  
}

# format
mod_dat <- mod_dat %>%
  rename(cluster = v001, years = time, urban = v025, Y = died)

# save
saveRDS(mod_dat, file = paste0("Data/", countries, "/nmr_data_prepped.rds"))


