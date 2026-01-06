# ============================================================================================================
#                                               --- START-UP ---
# ============================================================================================================


# Install packages for spatial data mapping and handling
install.packages(c(
  "dplyr",
  "ggplot2",
  "sysfonts",
  "showtext",
  "terra", 
  "geodata", 
  "rnaturalearth", 
  "rnaturalearthdata",
  "dismo", 
  "tidyverse", 
  "rgbif", 
  "raster" 
))

# Load packages
library(dplyr)
library(ggplot2)
library(terra)
library(geodata)
library(rnaturalearth)
library(dismo)
library(tidyr)
library(rgbif)
library(sysfonts)
library(showtext)

# Set working directory
wd <- setwd("C:/University Work/Year 3 Biology/HT/R Coding Assignment/Geo-data and Methods in R - Assigment")

# Create folders 
folders <- c("data/raw", "data/processed", "outputs/maps")
for (f in folders) {
  dir.create(file.path(wd, f), recursive = TRUE, showWarnings = FALSE)
}


# ============================================================================================================
#                                              --- TASK 1 ---
# ============================================================================================================

# This section involves running GLMs to  predict the present-day distribution of Species 1
# and Species 2 using climate variables. Then, these are used to present maps of their
# current distribution.
# Here, a function is written for the user to input the two species they are interested in,
# and what region they want to focus on (Europe, North/South America, Asia, Oceania, Africa).


# Region presets with 
REGION_PRESETS <- list(
  "Africa"        = c(-20, 55, -40, 40),
  "Europe"        = c(-25, 60, 34, 72),
  "South America" = c(-85, -30, -60, 15),
  "North America" = c(-170, -50, 5, 85),
  "Oceania"       = c(110, 180, -50, 0),
  "Asia"          = c(25, 180, -10, 85),
  "Global"        = c(-180, 180, -90, 90)
)


# Function generates current SDMs for a pair of species in a user-defined region
run_current_sdm <- function(species1, species2, region){
  
  # Create list to store the map outputs
  current_sdm_maps <- list()
  
  # Iterate current SDM generation for the user's species 1 and species 2
  chosen_species <- c(species1, species2)
  
  for (sp_name in chosen_species) {
    message(paste("Running current SDM generation for: ", sp_name))
    
    # Replace space " " with an underscore "_"
    sp_filename <- gsub(" ", "_", sp_name)
    
    
    # -1- Define file path and download occurrence data ------------------------------------------------------
    
    sp_file <- file.path("data", "raw", paste0(sp_filename, ".rds"))
    
    if (!file.exists(sp_file)) {
      message("Downloading occurrences from GBIF (may take a few seconds)...")
      # Limit download to 10,000 data points
      occ <- occ_search(scientificName = sp_name, hasCoordinate = TRUE, limit = 10000)
      occ_df <- occ$data
      saveRDS(occ_df, sp_file)
    } else {
      # Load if already exists
      occ_df <- readRDS(sp_file)
    }
    
    
    # -1.1- Extract and clean coordinates --------------------------------------------------------------------
    
    coords <- occ_df %>%
      dplyr::select(decimalLongitude, decimalLatitude) %>%
      rename(lon = decimalLongitude, lat = decimalLatitude)
    
    coords <- na.omit(coords)
    cat("Records with coordinates:", nrow(coords), "\n")
    
    
    
  }
}
  
  





































































