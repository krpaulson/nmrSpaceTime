
# Plot results

library(tidyverse)
library(magrittr)
library(SUMMER)
library(survey)
library(RColorBrewer)

# https://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
options(survey.lonely.psu="remove")


# settings ----------------------------------------------------------------

settings <- yaml::read_yaml("run_settings.yaml")
list2env(settings, .GlobalEnv)

years <- c(year_start_estimation:year_end_estimation)

survey_meta <- readRDS(paste0("Data/", country, "/metadata.rds"))
hold_out_year_start <- max(as.numeric(survey_meta$survey_year)) - 3
hold_out_years <- c(hold_out_year_start:2020)


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


# validation --------------------------------------------------------------

# use survey package to get direct estimates by admin1 and year, with SE
direct <- list()
for (y in unique(mod.dat$survey_year)) {
  design <- svydesign(
    ids = ~cluster,
    weights = ~v005,
    strata = ~strata,
    data = mod.dat[mod.dat$survey_year == y,]
  )
  direct[[y]] <- svyby(~Y, ~admin1.name + years, design, svymean)
  direct[[y]]$survey_year <- y
}
direct <- do.call(rbind, direct)

# light prep
direct <- direct %>%
  mutate(year = as.numeric(as.character(years)),
         admin1_name = admin1.name)

# combine with indirect and get MSE
outputMSE <- direct %>%
  left_join(results, by = c("admin1_name", "year")) %>%
  filter(year %in% hold_out_years & se > 0) %>%
  group_by(time_model) %>%
  mutate(w = 1/se^2) %>%
  summarise(mse = sum(w*(Y*1000 - nmr*1000)^2)/sum(w)) %>%
  arrange(mse)

# combine with indirect and get mean error
outputME <- direct %>%
  left_join(results, by = c("admin1_name", "year")) %>%
  filter(year %in% hold_out_years & se > 0) %>%
  group_by(time_model) %>%
  mutate(w = 1/se^2) %>%
  summarise(me = sum(w*(nmr*1000 - Y*1000))/sum(w)) %>%
  arrange(me)

# direct variance via delta method
direct$V <- with(direct, (1/Y + 1/(1-Y))^2 * se^2)

# get combined variance on logit scale
resultsCPO <- results %>%
  left_join(direct, by = c("admin1_name", "year")) %>%
  filter(year %in% hold_out_years & !is.na(Y)) %>%
  mutate(VarTot = V + logit_nmr_sd^2)

# compute CPO
resultsCPO$CPO_i <- 
  with(resultsCPO,
       -log(dnorm(x = SUMMER::logit(Y),
                  mean = SUMMER::logit(nmr),
                  sd = sqrt(VarTot))))

# summarize
outputCPO <- resultsCPO %>%
  filter(!is.na(CPO_i) & is.finite(CPO_i)) %>%
  group_by(time_model) %>%
  summarise(CPO = mean(CPO_i)) %>%
  arrange(CPO)

# combine results
output <- outputCPO %>%
  left_join(outputMSE, by = "time_model") %>%
  left_join(outputME, by = "time_model") %>%
  arrange(CPO)

# save
saveRDS(output, paste0("Results/", country, "/validation_results.rds"))


# plot -------------------------------------------------------------------

results$time_model <- factor(
  results$time_model,
  levels = c("ar1", "rw1", "rw2", "ar1_iv", "rw1_iv", "rw2_iv")
)

pdf(paste0("Results/", country, "/timeplot_summary.pdf"), width = 11, height = 6)
gg <- ggplot(results, aes(x = year, y = nmr * 1000, color = admin1_name, lty = admin1_name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Year", y = "NMR", color = "Admin 1", lty = "Admin 1", title = paste0(country, " NMR")) +
  scale_linetype_manual(
    values = rep(c("solid", "dashed", "dotted", "dotdash", "longdash",
                   "twodash", "F1"), 10)[1:length(unique(results$admin1_name))]) +
  scale_x_continuous(breaks = seq(2000, 2020, 5), minor_breaks = seq(2000, 2020, 5)) +
  facet_wrap("time_model")
print(gg)
dev.off()

pdf(paste0("Results/", country, "/timeplot_by_model.pdf"),
    width = 11, height = length(unique(results$admin1_name)))
for (tm in sort(unique(results$time_model))) {
  gg <- results %>%
    filter(time_model == tm) %>%
    ggplot(aes(x = year, y = nmr * 1000)) +
    geom_ribbon(aes(ymin = nmr_lower * 1000, ymax = nmr_upper * 1000), alpha = 0.4) +
    geom_line() +
    geom_point(data = direct, aes(y = Y * 1000, color = as.factor(survey_year))) +
    theme_bw() +
    labs(x = "Year", y = "NMR", title = paste0(country, " NMR // ", tm), color = "Survey Year") +
    scale_x_continuous(breaks = seq(2000, 2020, 5), minor_breaks = seq(2000, 2020, 5)) +
    facet_wrap("admin1_name", ncol = 3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(gg)
}
dev.off()

pdf(paste0("Results/", country, "/timeplot_by_admin1.pdf"), width = 11, height = 7)
for (adm1 in sort(unique(results$admin1_name))) {
  gg <- results %>%
    filter(admin1_name == adm1) %>%
    ggplot(aes(x = year, y = nmr * 1000)) +
    geom_ribbon(aes(ymin = nmr_lower * 1000, ymax = nmr_upper * 1000), alpha = 0.4) +
    geom_line() +
    geom_point(data = direct[direct$admin1.name == adm1,],
               aes(y = Y * 1000, color = as.factor(survey_year))) +
    theme_bw() +
    labs(x = "Year", y = "NMR", title = paste0(country, " NMR // ", adm1), color = "Survey Year") +
    scale_x_continuous(breaks = seq(2000, 2020, 5), minor_breaks = seq(2000, 2020, 5)) +
    facet_wrap("time_model", ncol = 3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(gg)
}
dev.off()



