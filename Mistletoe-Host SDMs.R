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
#                                           --- USER SET-UP ---
# ============================================================================================================

# Define species 1 and species 2
sp1 <- "Viscum album"
sp2 <- "Malus domestica"

# Define chosen bioclimatic variables for each species
sp_predictors <- list()
sp_predictors[[sp1]] <- c("bio1", "bio6", "bio11")
sp_predictors[[sp2]] <- c("bio6", "bio11", "bio10")

# Region presets with lon and lat bounds - add more if desired
REGION_PRESETS <- list(
  "Africa"        = c(-20, 55, -40, 40),
  "Europe"        = c(-25, 60, 34, 72),
  "South America" = c(-85, -30, -60, 15),
  "North America" = c(-170, -50, 5, 85),
  "Oceania"       = c(110, 180, -50, 0),
  "Asia"          = c(25, 180, -10, 85),
  "Global"        = c(-180, 180, -90, 90)
)


# ============================================================================================================
#                                              --- TASK 1 ---
# ============================================================================================================

# This section involves running GLMs to  predict the present-day distribution of Species 1
# and Species 2 using climate variables. Then, these are used to present maps of their
# current distribution.
# Here, a function is written for the user to input the two species they are interested in,
# and what region they want to focus on (Europe, North/South America, Asia, Oceania, Africa).

# Function generates current SDMs for a pair of species in a user-defined region
run_current_sdm <- function(species1, species2, region, predictor_list) {
  
  # Create list to store the map outputs
  current_sdm_maps <- list()
  
  # Iterate current SDM generation for the user's species 1 and species 2
  chosen_species <- c(species1, species2)
  
  for (sp_name in chosen_species) {
    message(paste("Running current SDM generation for: ", sp_name))
    
    # Replace space " " with an underscore "_"
    sp_filename <- gsub(" ", "_", sp_name)
    
    
    # -1- Define file path and download occurrence data ------------------------------------------------------
    
    sp_file <- file.path(wd, "data", "raw", paste0(sp_filename, ".rds"))
    
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
    
    
    # -1.2- Clip to study region chosen by user --------------------------------------------------------------
    
    if (!region %in% names(REGION_PRESETS)) {
      stop("Region not found. Please choose from the preset list.")
    }
    
    # Get boundaries for chosen region
    bounds <- REGION_PRESETS[[region]]
    coords_region <- coords %>%
      filter(lon >= bounds[1], lon <= bounds[2],
             lat >= bounds[3], lat <= bounds[4])
    
    # Create extent for chosen region
    region_ext <- ext(bounds[1], bounds[2], bounds[3], bounds[4])
    
    
    # -1.3- Remove occurrence points in the ocean ------------------------------------------------------------
    
    # Download the ocean data
    ocean_data_dir <- file.path(wd, "data", "raw", "ocean")
    
    if (!dir.exists(ocean_data_dir)) dir.create(ocean_data_dir)
    URL <- "https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_ocean.zip"
    zip_file <- file.path(ocean_data_dir, basename(URL))
    if (!file.exists(zip_file)) {
      download.file(URL, zip_file)
    }
    
    files <- unzip(zip_file, exdir = ocean_data_dir)
    
    # Find the shapefile (.shp)
    shp_file <- files[grepl("\\.shp$", files)]
    
    # Read with terra
    ocean <- vect(shp_file)
    
    # Convert coordinates to SpatVector
    species_vect <- vect(coords_region, geom = c("lon", "lat"), crs = "EPSG:4326")
    
    # Make sure CRS matches ocean
    crs(species_vect) <- crs(ocean)
    
    # Generate a matrix that checks every point against every ocean polygon
    ocean_intersects <- relate(species_vect, ocean, relation = "intersects")
    
    # Create a logical vector that returns TRUE if a point is in any ocean polygon
    is_ocean <- apply(ocean_intersects, 1, any)
    
    # Store the points that DO NOT intersect with the ocean polygons
    species_land_vect <- species_vect[!is_ocean, ]
    
    # Convert back to data.frame of coordinates
    species.coords <- as.data.frame(geom(species_land_vect)[, c("x", "y")])
    colnames(species.coords) <- c("lon", "lat")
    
    
    # -2- Environmental data and WorldClim bioclimatic variables ---------------------------------------------
    
    # Create file path for the climate data
    bioclim_dir <- file.path(wd, "data", "raw", "worldclim")
    if (!dir.exists(bioclim_dir)) dir.create(bioclim_dir)
    
    # Download 19 bioclim variables at 10 arc-minute resolution
    bioclim <- worldclim_global(var = "bio", res = 10, path = bioclim_dir)
    names(bioclim) <- paste0("bio", 1:19)
    
    
    # -2.1- Define a study extent with a buffer --------------------------------------------------------------
    
    xy <- crds(species_land_vect)
    study_ext <- ext(min(xy[,1]) - 5, max(xy[,1]) + 5, min(xy[,2]) - 5, max(xy[,2]) + 5)
    bioclim_crop <- crop(bioclim, study_ext)
    
    
    # -2.2- Extract raster values at occurrence points -------------------------------------------------------
    
    # Ensure point CRS matches rasters
    species_pts <- project(species_land_vect, crs(bioclim_crop))
    
    # Extract climate values at points
    clim_vals <- terra::extract(bioclim_crop, species_pts)
    
    # Combine coordinates with climate values
    species_data <- bind_cols(
      as_tibble(crds(species_pts)) %>% rename(lon = x, lat = y),
      as_tibble(clim_vals)[,-1]
    ) %>% drop_na()
    
    # Save processed data
    write.csv(species_data, file.path(wd,"data", "processed", paste0(sp_filename, ".csv")))
    
    
    # -3- Current Species Distribution Model Generation ------------------------------------------------------
    
    # -3.1- Generate background/pseudo-absence points -------------------------------------------------------------------
    
    # Sample random points from the first layer of the cropped climate data
    set.seed(123)
    bg_n <- 500
    bg_pts <- spatSample(bioclim_crop[[1]], size = bg_n, method = "random",
                         na.rm = TRUE, as.points = TRUE, values = FALSE)
    
    # Extract climate values for these background points
    bg_clim <- terra::extract(bioclim_crop, bg_pts)
    
    # Format background point data
    bg_coords <- as.data.frame(crds(bg_pts))
    colnames(bg_coords) <- c("lon", "lat")
    
    background_data <- bind_cols(bg_coords, as_tibble(bg_clim)[,-1]) %>%
      mutate(presence = 0) %>%
      drop_na()
    
    # Format presence data (Presence = 1)
    presence_data <- species_data %>% mutate(presence = 1)
    
    
    # -3.2- Split data into training (70%) and testing (30%) data
    
    # Split presence data
    k = 0.7
    pres_idx <- sample(nrow(presence_data), size = floor(k * nrow(presence_data)))
    train_pres <- presence_data[pres_idx, ]
    test_pres <- presence_data[-pres_idx, ]
    
    # Split background data
    bg_idx <- sample(nrow(background_data), size = floor(k * nrow(background_data)))
    train_bg <- background_data[bg_idx, ]
    test_bg <- background_data[-bg_idx, ]
    
    # Combine for training
    train_data <- bind_rows(train_pres, train_bg)
    
    
    # -3.3- Fit GLM -----------------------------------------------------------------------------------------------------
    
    # Identify variables for this species defined by the user
    current_vars <- predictor_list[[sp_name]]
    if(is.null(current_vars)) stop(paste("No variables stored for", sp_name))
    
    message(paste("Fitting GLM with: ", paste(current_vars, collapse = ", ")))
    
    # Construct formula dynamically
    formula_str <- paste("presence ~ ", paste(current_vars, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    # Fit model using training data
    sdm_model <- glm(model_formula, data = train_data, family = binomial)
    
    
    # -3.4- Evaluate model ----------------------------------------------------------------------------------------------
    
    # Evaluate using held-out test data
    eval_res <- evaluate(
      p = test_pres %>% dplyr::select(current_vars),
      a = test_bg %>% dplyr::select(current_vars),
      model = sdm_model
    )
    cat("AUC score:", eval_res@auc, "\n")
    
    
    # -3.5- Predict and Map ---------------------------------------------------------------------------------------------
    
    # Ensure raster has those layers
    if (!all(current_vars %in% names(bioclim_crop))) {
      stop("ERROR: One or more predictor layers not found in bioclim_crop: ", paste(pred_vars, collapse = ", "))
    }
    
    # Predict suitability across the full study extent
    prediction <- terra::predict(bioclim_crop[[current_vars]], sdm_model, type = "response")
    
    # Store in list to contain map using the species name as a label
    current_sdm_maps[[sp_name]] <- prediction
    
    # Plot
    plot(prediction, main = paste("SDM:", sp_name, "\nAUC:", round(auc_score, 2)))
    points(species_land_vect, pch = 16, cex = 0.5, col = "black")
    
    message("Completed generation of current SDM for:", sp_name)
    
    cat("--------------------------------------------------------------------------\n")

  }
  
  message("All SDMs generated successfully. Returning map list...")
  return(current_sdm_maps)
  
}

run_current_sdm(sp1, sp2, "Europe", sp_predictors)  
  





































































