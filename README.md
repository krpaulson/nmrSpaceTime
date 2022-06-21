# nmrSpaceTime

Estimate NMR by space and time, using DHS household survey data.

# Getting setup

1.  Register for a project to get access to DHS data: <https://dhsprogram.com/>

2.  Install the `rdhs` R package and configure your DHS project settings in R:

```{r, eval=F}
rdhs::set_rdhs_config(
  email = {YOUR_EMAIL},
  project = {YOUR_PROJECT_NAME}
)
```

3.  Clone this repository onto your computer.

4.  Open `nmrSpaceTime.Rproj` R project.

5.  Select a country and fill out your settings in `run_settings.yaml`. The settings are:

-   `country`: full name of the country (e.g. Zambia)
-   `country_code`: three-letter code for the country matching the country code in the DHS zip files. These are ISO3 codes which are standardized. (e.g. ZMB)
-   `gadmFolder`: folder name for GADM shape files (e.g. gadm40_ZMB; can find here: <https://gadm.org/download_country.html>)
-   `year_start_survey`: earliest DHS survey year using the most recent sampling frame
-   `frame_year`: sampling frame year, usually the last census, can find by searching "sampling frame" in the DHS report.

6.  Create the folder `nmrSpaceTime/Data/{country}`

7.  Create Admin 1 urban fraction file (`frame_urb_prop.csv` with columns "region" and "urban_prop"; following details here: <https://github.com/wu-thomas/SUMMER-DHS/tree/main/Vignette/Urban%20Fraction>) and put it in the country data folder.

# Running the code

Run the following scripts in order. No specific modifications should be required.

-   `0_prep.R`
-   `0b_prep_UR.R`
-   `1_fit.R`
-   `2_plot.R`
