# ============================================================================================
#                                       --- START-UP ---
# ============================================================================================

# Install packages for spatial data mapping and handling
install.packages(c(
  "dplyr",
  "ggplot2",
  "sysfonts", 
  "showtext",
  "here", 
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
library(here)
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
wd <- setwd("C:/University Work/Year 3 Biology/HT/R Coding Assignment/Geo-data and Methods in R - Assignment")

# Create folders 
folders <- c("data/raw", "data/processed", "outputs/maps")
for (f in folders) {
  dir.create(file.path(wd, f), recursive = TRUE, showWarnings = FALSE)
}


