
# Plot results

library(tidyverse)
library(magrittr)
library(SUMMER)


# settings ----------------------------------------------------------------

country <- "Liberia"
country_code <- "LBR"
years <- c(2000:2020)

hold_out_years <- c(2009:2011)


# setup -------------------------------------------------------------------

# load GADM
dsn <- paste0("Data/", country, "/shapeFiles_gadm")
layer <- list.files(dsn)
layer <- layer[grep("_1.shp", layer)]
layer <- gsub(".shp", "", layer)
poly_adm1 <- rgdal::readOGR(dsn = dsn, layer = layer)

# load data
mod.dat <- readRDS(paste0("Data/", country, "/nmr_data_prepped.rds"))

# load results

results <- list.files(
  paste0("Results/", country),
  pattern = "results_holdout_",
  full.names = T
) %>% map_dfr(readRDS)

results_all <- readRDS(paste0("Results/", country, "/results_all.rds"))


# validation ---------------------------------------------------------------

# direct
direct <- mod.dat %>%
  mutate(year = as.numeric(as.character(years)),
         admin1_name = admin1.name) %>%
  group_by(admin1_name, admin1.char, year, survey_year) %>%
  summarise(direct = sum(Y*v005)/sum(v005)) %>%
  ungroup()

# combine with indirect and get MSE
direct %>%
  left_join(results, by = c("admin1_name", "year")) %>%
  filter(year %in% hold_out_years) %>%
  group_by(time_model) %>%
  summarise(mse = mean((direct - nmr)^2))

# get CPO (log score)


# plot -------------------------------------------------------------------

pdf(paste0("Results/", country, "/timeplot_summary.pdf"), width = 11, height = 6)
ggplot(results_all, aes(x = year, y = nmr * 1000, color = admin1_name, lty = admin1_name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Year", y = "NMR", color = "Admin 1", lty = "Admin 1", title = "Liberia NMR") +
  scale_linetype_manual(
    values = rep(c("solid", "dashed", "dotted", "dotdash", "longdash",
                   "twodash", "F1"), 3)[1:length(unique(results_all$admin1_name))]) +
  scale_x_continuous(breaks = seq(1990, 2025, 5), minor_breaks = seq(1990, 2025, 5)) +
  facet_wrap("time_model")
dev.off()

pdf(paste0("Results/", country, "/timeplot_by_model.pdf"), width = 11, height = 10)
for (tm in sort(unique(results_all$time_model))) {
  gg <- results_all %>%
    filter(time_model == tm) %>%
    ggplot(aes(x = year, y = nmr * 1000)) +
    geom_ribbon(aes(ymin = nmr_lower * 1000, ymax = nmr_upper * 1000), alpha = 0.4) +
    geom_line() +
    geom_point(data = direct, aes(y = direct * 1000, color = as.factor(survey_year))) +
    theme_bw() +
    labs(x = "Year", y = "NMR", title = paste0("Liberia NMR // ", tm), color = "Survey Year") +
    scale_x_continuous(breaks = seq(1990, 2020, 5), minor_breaks = seq(1990, 2020, 5)) +
    scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 25)) +
    facet_wrap("admin1_name", ncol = 3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(gg)
}
dev.off()

pdf(paste0("Results/", country, "/CIs-2019.pdf"))
adm1_order <- direct %>%
  filter(time_model == "rw2_iv" & year == 2019) %>%
  arrange(nmr) %>%
  dplyr::select(admin1_name)
for (tm in sort(unique(results_all$time_model))) {
  gg <- results_all %>%
    filter(time_model == tm & year == 2019) %>%
    arrange(nmr) %>%
    mutate(admin1_name = factor(admin1_name, levels=adm1_order$admin1_name)) %>%
    ggplot(aes(y = admin1_name)) +
    geom_errorbarh(aes(xmin = nmr_lower * 1000, xmax=nmr_upper * 1000), height = 0.4) +
    geom_point(aes(x = nmr * 1000)) +
    geom_point(data = filter(direct, year >= 2015), aes(x = direct * 1000, y = admin1_name), color = "magenta") +
    theme_bw() +
    labs(x = "NMR", y = "", title = paste0("Liberia NMR 2019 // ", tm)) +
    scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20), minor_breaks = seq(0, 80, 20))
  print(gg)
}
dev.off()


# Maps --------------------------------------------------------------------

expit_medians$nmr1000 <- expit_medians$nmr * 1000
expit_medians <- as.data.frame(expit_medians)

expit_medians_subset <- 
  expit_medians[expit_medians$year %in% seq(2000, 2020, 5),]

pdf(paste0("Results/", country, "/map_adm1.pdf"))
SUMMER::mapPlot(data = expit_medians_subset,
                geo = poly_adm1,
                variables = "year",
                values = "nmr1000",
                by.data = "admin1_name",
                by.geo = "NAME_1",
                legend.label = "NMR",
                is.long = T)
dev.off()


