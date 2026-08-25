# ==============================================================================
# Confirm species presence 
#
# Description: Shared first step for this whole set of scripts. Takes a raw
# acoustic classifier output and a manual verification dataset, and returns
# the raw detections restricted to site-species pairs with at least one
# confirmed TRUE detection (matched by exact file name).
#
# Required input columns:
#   raw data:          file, species, score
#   verification data: file, species, identity (TRUE/FALSE)
#
# NOTE: site (recorder location) is extracted by selecting everything before
# the first underscore in the file name - change this to match your string.
# ==============================================================================
 
required_packages <- c("dplyr", "stringr", "readr")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)
 
library(dplyr)
library(stringr)
library(readr)
 
get_confirmed_detections <- function(data, verification) {
 
  df <- data %>% mutate(site = str_extract(file, "^[^_]+"))
 
  confirmed <- verification %>%
    filter(identity == TRUE | identity == "TRUE" | identity == "['TRUE']") %>%
    distinct(file, species)
 
  confirmed_pairs <- df %>%
    inner_join(confirmed, by = c("file", "species")) %>%
    distinct(site, species)
 
  df %>% semi_join(confirmed_pairs, by = c("site", "species"))
}
 
# ==============================================================================
# Example usage
# ==============================================================================
 
# raw <- read_csv("data/classifier_output.csv")
# ver <- read_csv("data/verification_file.csv")
#
# confirmed <- get_confirmed_detections(raw, ver)
