# ==============================================================================
# Calculate daily acoustic activity from a classifier output
#
# Description: Takes a single acoustic classifier output (e.g. BTO, BirdNet,
# Polyglotta) and calculates daily activity (number of detections) per species 
# per site (recorder location). Function used in phenology.R and
# restoration_comparison.R.
#
# Requires: source("code/verify_presence.R") first, and a `confirmed` object
# from get_confirmed_detections(raw, ver). calculate_daily_activity() takes
# BOTH the original raw data AND confirmed - not confirmed alone - because
# working out which site-days were actually surveyed (step 1 below) needs
# the FULL raw data. If it only saw the confirmed subset, a site-day where
# the recorder ran but no CONFIRMED species happened to call that day would
# disappear entirely, breaking the true-zero-fill this function is built
# around.
#
# Required input columns:
#   data:      file, species, score
#   confirmed: file, species, site (i.e. the output of get_confirmed_detections())
#
# NOTE: site (recorder location) is extracted by selecting everything before
# the first underscore in the file name - change this to match your string.
#
# Confidence filtering is controlled by a single `threshold` argument, which
# accepts either:
#   - one number, applied to every species, or
#   - a data frame with columns `species` and `threshold`, giving each
#     species its own cutoff. Any species NOT listed in that data frame gets
#     no cutoff and is dropped from the results entirely (it will show 0
#     activity everywhere). See species_specific_thresholds.R for how to
#     derive these thresholds from validated data.
# ==============================================================================

# Install any missing packages, then load them
required_packages <- c("dplyr", "lubridate", "stringr", "tidyr", "readr")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(readr)

# Function to calculate daily activity per species per site
calculate_daily_activity <- function(data,
                                     confirmed,
                                     tz,
                                     threshold) {
  
  # --- 1. Work out which site-days were actually surveyed ---
  # This must come from the FULL raw dataset, before any species filtering,
  # since survey effort (was a recorder running that day) is independent
  # of which species turn out to be confirmed.
  df <- data %>%
    mutate(
      site = str_extract(file, "^[^_]+"),
      datetime = ymd_hms(str_extract(file, "\\d{8}[_-]\\d{6}"), tz = tz),
      date = as_date(datetime)
    )
  
  survey_days <- df %>% distinct(site, date)
  
  # --- 2. Restrict to confirmed site-species pairs ---
  confirmed_pairs <- confirmed %>% distinct(site, species)
  df <- df %>% semi_join(confirmed_pairs, by = c("site", "species"))
  
  all_species <- unique(confirmed$species)
  
  # --- 3. Apply confidence threshold(s) ---
  if (is.data.frame(threshold)) {
    # Species-specific: join each species to its own cutoff. Using
    # inner_join (not left_join) means any species missing from the
    # threshold table is dropped here rather than silently falling back to
    # some default - see the note at the top of this file.
    df_filtered <- df %>%
      inner_join(threshold, by = "species") %>%
      filter(score >= threshold)
  } else if (is.numeric(threshold) && length(threshold) == 1) {
    df_filtered <- df %>% filter(score >= threshold)
  } else {
    stop("`threshold` must be a single number, or a data frame with columns `species` and `threshold`")
  }
  
  # --- 4. Summarise activity and fill in explicit zeros ---
  daily_activity <- survey_days %>%
    expand_grid(species = all_species) %>%
    left_join(
      df_filtered %>% count(site, date, species, name = "activity"),
      by = c("site", "date", "species")
    ) %>%
    mutate(activity = replace_na(activity, 0)) %>%
    arrange(site, date, species)
  
  return(daily_activity)
}

# ==============================================================================
# Example usage
# ==============================================================================

# source("code/verify_presence.R")
#
# raw <- read_csv("data/classifier_output.csv")
# ver <- read_csv("data/verification_file.csv")
# confirmed <- get_confirmed_detections(raw, ver)
#
# # Option 1: one threshold applied to every species
# daily_activity <- calculate_daily_activity(raw, confirmed, tz = "UTC", threshold = 0.5)
#
# # Option 2: species-specific thresholds - any species NOT listed here is
# # dropped from the results (shows 0 activity everywhere)
# source("code/species_specific_thresholds.R")
# sp_thresholds <- calculate_precision_thresholds(ver, 0.9)$thresholds
# daily_activity <- calculate_daily_activity(raw, confirmed, tz = "UTC",
#                                            threshold = sp_thresholds)
