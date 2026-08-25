# ==============================================================================
# Species-specific confidence thresholds from verified detections
#
# Description: Takes a set of manually verified detections (each with a
# confidence score and a TRUE/FALSE identity) and, for each species, fits a
# logistic curve of score against precision (probability a detection is
# correct). Returns the minimum score needed to reach each target precision
# (e.g. 0.5, 0.7, 0.9).
#
# Required input columns: file, species, identity (TRUE/FALSE)
#
# NOTE: this needs a reasonable number of validated detections per species
# (spread across the confidence range) to fit a sensible curve - a handful
# of points, or points all clustered at one end of the score range, will
# give an unstable or meaningless threshold.
# ==============================================================================

# Install any missing packages, then load them
required_packages <- c("dplyr", "purrr", "tidyr", "ggplot2", "readr")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(readr)

# Function to calculate species-specific confidence scores for target precisions
calculate_precision_thresholds <- function(data,
                                           precisions) {
  
  # --- 1. Create a binary truth variable ---
  df <- data %>%
    mutate(truth = if_else(
      identity == TRUE | identity == "TRUE" | identity == "['TRUE']", 1, 0
    ))
  
  # --- 2. Fit a logistic curve of score -> precision, per species ---
  species_list <- unique(df$species)
  
  curves <- map_dfr(species_list, function(sp) {
    sp_data <- df %>% filter(species == sp)
    
    # Skip species with too little data or no variation in truth (glm will
    # fail or give a meaningless curve)
    if (nrow(sp_data) < 10 || n_distinct(sp_data$truth) < 2) {
      warning("Skipping ", sp, ": not enough validated data to fit a curve")
      return(NULL)
    }
    
    mod <- glm(truth ~ score, data = sp_data, family = binomial)
    
    tibble(
      species = sp,
      score = seq(0, 1, by = 0.01),
      precision = predict(mod, newdata = tibble(score = seq(0, 1, by = 0.01)),
                          type = "response")
    )
  })
  
  # --- 3. For each species and target precision, find the minimum score
  #        at which the fitted curve first reaches that precision ---
  thresholds <- map_dfr(precisions, function(p) {
    curves %>%
      group_by(species) %>%
      filter(precision >= p) %>%
      slice_min(score, n = 1) %>%
      ungroup() %>%
      transmute(species, precision = p, threshold = score)
  })
  
  list(thresholds = thresholds, curves = curves)
}

# ==============================================================================
# Example usage
# ==============================================================================

# ver <- read_csv("data/verification_file.csv")
#
# result <- calculate_precision_thresholds(ver, c(0.5, 0.7, 0.9))
# result$thresholds
#
# # Bin raw detections and calculate mean number of TRUE detections
# points <- ver %>%
#   mutate(truth = if_else(
#     identity == TRUE | identity == "TRUE" | identity == "['TRUE']", 1, 0
#   ),
#   score_bin = round(score, 1)) %>%
#   group_by(species, score_bin) %>%
#   summarise(mean_truth = mean(truth), n = n(), .groups = "drop")
#
# # Plot logistic curves (one panel per species)
# ggplot() +
#     geom_point(data = points,
#                aes(x = score_bin, y = mean_truth, size = n),
#                alpha = 0.4) +
#     geom_line(data = result$curves, aes(x = score, y = precision), linewidth = 1) +
#     geom_vline(data = result$thresholds %>% filter(precision == 0.9),
#                aes(xintercept = threshold), linetype = "dashed", colour = "red") +
#     facet_wrap(vars(species)) +
#     labs(x = "Confidence score", y = "Precision", size = "# checked") +
#     theme_classic()
