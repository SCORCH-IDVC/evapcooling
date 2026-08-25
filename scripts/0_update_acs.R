library(here)
library(sf)
library(tidycensus)

# Back up the original (2017-2021 ACS)
original <- here("data", "Q1 Data Shapefile", "pima_Q1_data.shp")
backup <- here("data", "Q1 Data Shapefile", "pima_Q1_data_acs2017_2021.shp")

if (!file.exists(backup)) {
  bg_orig <- st_read(original)
  st_write(bg_orig, backup)
  cat("Original backed up to:", backup, "\n")
} else {
  cat("Backup already exists. Skipping.\n")
}

# Read the current shapefile
bg_sf <- st_read(original)

# Pull 2019-2023 ACS5 for Pima County block groups
acs <- get_acs(
  geography = "block group",
  state = "AZ",
  county = "Pima",
  year = 2023,
  survey = "acs5",
  variables = c(
    med_income = "B19013_001",
    pop_total = "B03002_001",
    pop_white_nh = "B03002_003"
  ),
  output = "wide"
)

acs$pct_minority <- 1 - (acs$pop_white_nhE / acs$pop_totalE)

acs_clean <- data.frame(
  GEOID = sub("^0", "", acs$GEOID),
  med_inc_new = acs$med_incomeE,
  pct_mnr_new = acs$pct_minority,
  pop_new = acs$pop_totalE
)

# 4. Match and replace
bg_sf$geoid20 <- sub("^0", "", bg_sf$geoid20)
bg_sf <- merge(bg_sf, acs_clean, by.x = "geoid20", by.y = "GEOID", all.x = TRUE)

bg_sf$med_inc <- bg_sf$med_inc_new
bg_sf$pct_mnr <- bg_sf$pct_mnr_new
bg_sf$pop <- bg_sf$pop_new
bg_sf$med_inc_new <- NULL
bg_sf$pct_mnr_new <- NULL
bg_sf$pop_new <- NULL

# Overwrite the default shapefile
st_write(bg_sf, original, delete_dsn = TRUE)
cat("Updated shapefile written to:", original, "\n")
cat("All scripts will now use 2019-2023 ACS data.\n")
