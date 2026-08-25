# ==============================================================================
# Exploring phenology from processed acoustic activity data
#
# Description: Takes daily activity data (one row per site/date/species with an
# `activity` count) and:
#   1. Smooths each species' activity with a rolling mean over a chosen
#      window (to remove day-to-day noise), then scales it 0-1 (so loud
#      and quiet species become directly comparable).
#   2. Plots a heatmap of scaled activity over time, one row per species,
#      ordered by the date each species' activity peaks - a quick way to
#      see phenology (which species are active when) across a season.
#
# Required input columns:
#   raw data:          file, species, score
#   verification data: file, species, identity (TRUE/FALSE)
#
# NOTE: this pools activity across all sites - site isn't part of the
# grouping anywhere in this script, so a species' smoothed value on a given
# date reflects activity summed/averaged across every site that recorded
# that day. To adapt this for a by-site breakdown instead, look for the
# three lines marked "<- for by-site" in generate_smoothed_dataset() below
# (two group_by() calls plus the final select()) - all three need `site`
# added alongside `species`. You'd also want to carry `site` through to
# plot_phenology_heatmap() if you want it distinguished on the plot too
# (e.g. via facet_wrap(~ site)).
#
# NOTE: site (recorder location) is extracted by selecting everything before
# the first underscore in the file name - change this to match your string.
# ==============================================================================

# Install any missing packages, then load them
required_packages <- c("dplyr", "lubridate", "slider", "ggplot2", "readr")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(dplyr)
library(lubridate)
library(slider)
library(ggplot2)
library(readr)

source("code/verify_presence.R") # Load function to select only detections where
                                 # species has been manually verified as present
source("code/daily_activity_processing.R") # Load function to calculate daily 
                                           # activity

# Load datasets
raw <- read_csv("data/classifier_output.csv")
ver <- read_csv("data/verification_file.csv")
confirmed <- get_confirmed_detections(raw, ver)

# Calculate daily activity - one threshold applied to every species
daily_activity <- calculate_daily_activity(raw, confirmed, "UTC", 0.5)

# Select window to calculate rolling average
window <- 7
half_window <- floor(window / 2)
  
# Just in case, drop species with no activity at all - these would just create 
# flat, uninformative rows once scaled. Grouping by species here matters: it
# scopes the all() check to each species' own rows, so a species that's
# genuinely inactive everywhere gets dropped without touching any other
# species. Without the grouping, all() would check the WHOLE dataset at
# once and this filter would do essentially nothing.
daily_activity <- daily_activity %>%
  group_by(species) %>%  # <- for by-site: group_by(species, site)
  filter(!all(activity == 0 | is.na(activity))) %>%
  ungroup()

# Smooth, scale, and assemble the per-species phenology dataset
smoothed <-  daily_activity %>%
  group_by(species) %>%  # <- for by-site: group_by(species, site)
  arrange(date) %>%
  mutate(
    # Rolling mean, sliding over row INDICES within each group (not the
    # values directly) - this avoids instability that slide_index_dbl()
    # can show when applied straight to a column inside a grouped mutate.
    #
    # Each window also has a "gap check": if the dates actually present
    # in a window fall well short of the requested window width (e.g. a
    # recorder was down for a few days), the average is discarded rather
    # than calculated from just a couple of days - this avoids false
    # "peaks" appearing near the start/end of a series purely because
    # fewer days went into the average there.
    activity_smooth = {
      val_vec <- activity
      day_vec <- date
      
      slide_index_dbl(
        .x = seq_along(val_vec),
        .i = day_vec,
        .f = function(indices) {
          w_vals <- val_vec[indices]
          w_days <- day_vec[indices]
          span <- as.numeric(difftime(max(w_days), min(w_days), units = "days"))
          if (span < (window - 2)) return(NA_real_)  # allow one missing day
          mean(w_vals, na.rm = TRUE)
        },
        .before = days(half_window),
        .after  = days(half_window),
        .complete = TRUE
      )
    },
    
    # Min-max scale 0-1 so different species become comparable.
    # If every smoothed value in the group is NA (not enough days to
    # compute a single window), or the group is a flat line, handle those
    # cases explicitly rather than letting min()/max() error or return Inf.
    activity_scaled = if (all(is.na(activity_smooth))) {
      NA_real_
    } else {
      mn <- min(activity_smooth, na.rm = TRUE)
      mx <- max(activity_smooth, na.rm = TRUE)
      rng <- mx - mn
      if (rng == 0) 0 else (activity_smooth - mn) / rng
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(activity_scaled)) %>%
  select(species, date, activity_scaled) %>%  # <- for by-site: add site here too
  distinct()

# Plot a heatmap of scaled activity over time, one row per species, ordered by
# the date each species' activity peaks - a quick way to see phenology (which
# species are active when) across a season.
peak_order <- smoothed %>%
  group_by(species) %>%
  summarise(
    peak = date[which.max(activity_scaled)],
    total_activity = sum(activity_scaled, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(peak, desc(total_activity))

smoothed %>%
  mutate(species = factor(species, levels = peak_order$species)) %>%
  ggplot(aes(x = date, y = species)) +
  geom_tile(aes(fill = activity_scaled))
