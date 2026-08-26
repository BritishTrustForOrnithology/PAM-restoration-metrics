# ==============================================================================
# Acoustic metrics to monitor restoration progress
#
# Description: Calculates multiple ecological metrics from confirmed acoustic 
# detections, including:
#   1. Species lists (project-wide and per-site)
#   2. Species richness per site
#   3. Top-5-day activity means per site-species
#   4. Simpson's diversity index (1-D) per site
#   5. Occupancy model
#
# Required input columns:
#   raw data:          file, species, score
#   verification data: file, species, identity (TRUE/FALSE)
#
# NOTE: site (recorder location) is extracted by selecting everything before
# the first underscore in the file name - change this to match your string.
#
# Additional data input:
#   site metadata:     site, any covariates, lat/lon
# ==============================================================================

# Install any missing packages, then load them
required_packages <- c("dplyr", "lubridate", "stringr", "readr", "tidyr", "ggplot2", 
                       "tibble", "vegan", "sf", "spOccupanc")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(dplyr)
library(lubridate)
library(stringr)
library(readr)
library(tidyr)
library(ggplot2)
library(tibble)
library(vegan)
library(sf)
library(spOccupancy)

source("code/verify_presence.R") # Load function to select only detections where
                                 # species has been manually verified as present
source("code/daily_activity_processing.R") # Load function to calculate daily 
                                           # activity

# Load data 
raw <- read_csv("data/classifier_output.csv")
ver <- read_csv("data/verification_file.csv")
site_metadata <- read_csv("data/site_metadata.csv")
confirmed <- get_confirmed_detections(raw, ver)

# Project-wide species list ----

study_list <- confirmed %>%
  group_by(species) %>%
  summarise(
    detections = n(),
    sites = n_distinct(site),
    .groups = "drop"
  ) %>%
  select(species, sites, detections) %>%
  arrange(species)

# Per-site species list ----

site_list <- confirmed %>%
  count(site, species, name = "detections") %>%
  arrange(desc(detections))

# Richness per site ----

# NOTE: Ensure survey effort (e.g., number of days) is roughly equal across 
# sites. If highly unequal, consider rarefied richness instead.
richness_df <- confirmed %>%
  group_by(site) %>%
  summarise(richness = n_distinct(species), .groups = "drop")

# Mean activity of the top five days of greatest activity ----

# Calculate daily activity using a standard 0.5 threshold
# (Replace with species-specific thresholds if preferred)
daily_activity <- calculate_daily_activity(raw, confirmed, tz = "UTC", threshold = 0.5)

# Average over top 5 days
top5_activity <- daily_activity %>%
  group_by(site, species) %>%
  slice_max(activity, n = 5, with_ties = FALSE) %>%
  summarise(top5_mean = mean(activity), .groups = "drop") %>%
  # Filter out sites where the species was never genuinely confirmed present
  semi_join(confirmed %>% distinct(site, species), by = c("site", "species"))

# Example plot comparing activity between management types
top5_activity %>%
  left_join(site_management, by = "site") %>%
  ggplot(aes(x = management, y = top5_mean, fill = management)) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(vars(species), scales = "free_y") 

# Simpson's diversity index ----

# Build a site x species abundance matrix from TOTAL activity
site_matrix <- daily_activity %>%
  group_by(site, species) %>%
  summarise(total_activity = sum(activity), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = total_activity, values_fill = 0) %>%
  column_to_rownames("site")

# Drop empty sites to prevent vegan::diversity() from erroneously returning a 
# Simpson's index of 1 (the maximum possible value) for a site with 0 detections.
empty_sites <- rownames(site_matrix)[rowSums(site_matrix) == 0]
if (length(empty_sites) > 0) {
  warning(
    "Site(s) with zero total activity excluded from diversity calculation: ",
    paste(empty_sites, collapse = ", ")
  )
  site_matrix <- site_matrix[rowSums(site_matrix) > 0, , drop = FALSE]
}

# Calculate index
simpson_df <- tibble(
  site = rownames(site_matrix),
  simpson = diversity(site_matrix, index = "simpson")
)

# Occupancy model ---- 

# Get unique lists for array dimensions
sp_names <- sort(unique(confirmed$species))
unique_sites <- sort(unique(confirmed$site))

# Define repeated measures to calculate detectability
# Only keeps (site, date) pairs where recording actually took place
survey_effort <- raw %>%
  mutate(site = str_extract(file, "^[^_]+"),
         datetime = ymd_hms(str_extract(file, "\\d{8}[_-]\\d{6}"), tz = "UTC"),
         date = as_date(datetime)) %>%
  select(site, date) %>%
  distinct() %>%
  # Explicit chronological order before numbering occasions - row_number()
  # without this depends on incidental row order in `raw`, which usually
  # happens to be date-sorted but isn't guaranteed to be
  arrange(site, date) %>%
  group_by(site) %>%
  mutate(
    occasion = row_number(), # Assigns 1, 2, 3, etc. for each visit to a site
    julian = yday(date)      # Calculates Julian day (1-365) for your detection covariate
  ) %>%
  ungroup()

# Calculate number of repeats, species, and sites
n_occasions <- max(survey_effort$occasion)
n_sp <- length(sp_names)
n_sites <- length(unique_sites)

# Generate 1/0 grid for all survey days, joined to occasion number
daily_grid <- survey_effort %>%
  select(site, date, occasion) %>%
  # Cross active site-dates with all unique species
  expand_grid(species = sp_names) %>%
  left_join(
    # Distinct species presences per site per date (confirmed already has
    # `site` - no need to re-derive it from `file` again here)
    confirmed %>%
      mutate(datetime = ymd_hms(str_extract(file, "\\d{8}[_-]\\d{6}"), tz = "UTC"),
             date = as_date(datetime)) %>%
      select(site, date, species) %>%
      distinct() %>%
      mutate(detected = 1),
    by = c("site", "date", "species")
  ) %>%
  # Set 0 for species missing on an active survey date
  mutate(detected = replace_na(detected, 0))

# Initialize arrays with NA (NA means the site was not surveyed on that occasion)
y <- array(NA, dim = c(n_sp, n_sites, n_occasions))
julian_day <- array(NA, dim = c(n_sites, n_occasions))

# Name the dimensions for clarity
dimnames(y) <- list(
  species = sp_names,
  site = unique_sites,
  occasion = 1:n_occasions
)

# Fill the arrays with vectorized matrix-indexing, instead of looping over
# every (site, occasion, species) combination and re-scanning `confirmed`
# from scratch for each one (a filter() %>% nrow() > 0 check inside a triple
# loop does a full table scan per cell - fine for a handful of sites, but
# scales very badly, and daily_grid above already computed exactly this 0/1
# grid in one vectorized join, so there's no need to recompute it cell by
# cell).
julian_day[cbind(match(survey_effort$site, unique_sites), survey_effort$occasion)] <-
  survey_effort$julian

y[cbind(
  match(daily_grid$species, sp_names),
  match(daily_grid$site, unique_sites),
  daily_grid$occasion
)] <- daily_grid$detected

# Prepare covariates and run model
# Occupancy covariates (site-level)
occ_covs <- site_metadata %>%
  arrange(site) %>%
  select(site, forest_type, age)

# Detection covariates (observation-level)
det_covs <- list(
  julian = julian_day,
  site_id = matrix(1:n_sites, nrow = n_sites, ncol = n_occasions)
)

# Spatial coordinates
coords_proj <- site_metadata %>%
  arrange(site) %>%
  st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(2213) %>% # assuming project area is Caprathia, Romania
  st_coordinates()

# Bundle into spOccupancy format
sp_occ_df <- list(
  y = y,
  occ.covs = occ_covs,
  det.covs = det_covs,
  coords = coords_proj
)

# Define model parameters
# (Kept small here for a quick "simple" test run - often need long runs)
samples <- 10000
burn <- 5000
formula <- ~ factor(forest_type) + factor(age)
inits <- list(
  beta = 0, beta.comm = 0,
  tau.sq.beta = 0.5, sigma.sq.mu = 0.5, kappa = 0.5
)

# Run multi-species occupancy model
occ_out <- msPGOcc(
  occ.formula = formula,
  det.formula = ~ scale(julian) + scale(I(julian^2)) + (1 | site_id),
  data = sp_occ_df,
  inits = inits,
  n.samples = samples,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 1000,
  n.burn = burn,
  n.thin = 10,
  n.chains = 3
)

summary(occ_out)
