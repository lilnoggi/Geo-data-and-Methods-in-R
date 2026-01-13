# ==============================================================================================================
#                                               --- START-UP ---
# ==============================================================================================================


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
library(ggplot2)
library(here)
library(terra)
library(geodata)
library(rnaturalearth)
library(dismo)
library(tidyr)
library(rgbif)
library(sysfonts)
library(showtext)

# Create folders 
folders <- c("data/raw", "data/processed", "outputs/maps")
for (f in folders) {
  dir.create(here(f), recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================================================
#                                           --- USER SET-UP ---
# ==============================================================================================================

# Define species 1 and species 2
sp1 <- "Loranthus europaeus"
sp2 <- "Quercus petraea"

# Define chosen bioclimatic variables for each species
sp_predictors <- list()
sp_predictors[[sp1]] <- c("bio3", "bio8", "bio15")
sp_predictors[[sp2]] <- c("bio1", "bio4", "bio12")

# Region presets with lon and lat bounds - add more if desired
REGION_PRESETS <- list(
  "Africa"        = c(-20, 55, -40, 37.5),
  "Europe"        = c(-25, 60, 34, 72),
  "South America" = c(-95, -30, -60, 15),
  "North America" = c(-170, -50, 5, 85),
  "Oceania"       = c(110, 180, -50, 0),
  "Asia"          = c(25, 180, -10, 85),
  "Global"        = c(-180, 180, -90, 90)
)


# ==============================================================================================================
#                                              --- TASK 1 ---
# ==============================================================================================================

# This section involves running GLMs to  predict the present-day distribution of Species 1
# and Species 2 using climate variables. Then, these are used to present maps of their
# current distribution.
# Here, a function is written for the user to input the two species they are interested in,
# and what region they want to focus on (Europe, North/South America, Asia, Oceania, Africa).


# Function generates current SDMs for a pair of species in a user-defined region
run_current_sdm <- function(species1, species2, region, predictor_list) {
  
  # -1- Load static data once and prepare places to store outputs ----------------------------------------------
  
  message("Loading environmental and ocean data (may take a few seconds)...")
  
  
  # -1.1- Download the ocean data ------------------------------------------------------------------------------
  
  ocean_data_dir <- here("data", "raw", "ocean")
  
  if (!dir.exists(ocean_data_dir)) dir.create(ocean_data_dir)
  URL <- "https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_ocean.zip"
  zip_file <- file.path(ocean_data_dir, basename(URL))
  if (!file.exists(zip_file)) {
    download.file(URL, zip_file)
  }
  
  # Only unzips if it hasn't been unzipped already
  if (length(list.files(ocean_data_dir, pattern = "\\.shp$")) == 0) {
    unzip(zip_file, exdir = ocean_data_dir)
  }
  
  # Find the shapefile
  shp_file <- list.files(ocean_data_dir, pattern = "\\.shp$", full.names = TRUE)[1]
  
  # Read with terra
  ocean <- vect(shp_file)
  
  
  # -1.2- Download environmental data --------------------------------------------------------------------------
  
  # Create file path for the climate data
  bioclim_dir <- here("data", "raw", "worldclim")
  if(!dir.exists(bioclim_dir)) dir.create(bioclim_dir)
  
  # Download 19 bioclim variables at 10 arc-minute resolution
  bioclim_global <- worldclim_global(var = "bio", res = 10, path = bioclim_dir)
  names(bioclim_global) <- paste0("bio", 1:19)
  
  
  # -1.3- Prepare places to store outputs after for loop executed ----------------------------------------------
  
  # Create list to store the map outputs
  current_sdm_maps <- list()
  
  # Create a table to store results
  results_stats <- data.frame(
    Species = character(),
    Region = character(),
    N_Points = numeric(),
    AUC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate current SDM generation for the user's species 1 and species 2
  chosen_species <- c(species1, species2)
  
  for (sp_name in chosen_species) {
    message(paste("Running current SDM generation for: ", sp_name, " in ", region))
    
    # Replace space " " with an underscore "_"
    sp_filename <- gsub(" ", "_", sp_name)
    
    
    # -2- Define file path and download occurrence data --------------------------------------------------------
    
    sp_file <- here("data", "raw", paste0(sp_filename, ".rds"))
    
    if (!file.exists(sp_file)) {
      message("Downloading occurrences from GBIF (may take a few minutes)...")
      # Limit download to 10,000 data points
      occ <- occ_search(scientificName = sp_name, hasCoordinate = TRUE, limit = 10000)
      occ_df <- occ$data
      saveRDS(occ_df, sp_file)
    } else {
      # Load if already exists
      occ_df <- readRDS(sp_file)
    }
    
    # [SAFETY CHECK] Checks if there are any records at all
    if (is.null(occ_df) || nrow(occ_df) == 0) {
      message(paste("SKIP: No records found for", sp_name, "- check spelling or internet connection!"))
      next 
    }
    
    
    # -2.1- Extract and clean coordinates ----------------------------------------------------------------------
    
    coords <- occ_df %>%
      dplyr::select(decimalLongitude, decimalLatitude) %>%
      rename(lon = decimalLongitude, lat = decimalLatitude)
    
    coords <- na.omit(coords)
    cat("Records with coordinates:", nrow(coords), "\n")
    
    
    # -2.2- Clip to study region chosen by user ----------------------------------------------------------------
    
    if (!region %in% names(REGION_PRESETS)) {
      stop("Region not found. Please choose from the preset list.")
    }
    
    # Get boundaries for chosen region
    bounds <- REGION_PRESETS[[region]]
    coords_region <- coords %>%
      filter(lon >= bounds[1], lon <= bounds[2],
             lat >= bounds[3], lat <= bounds[4])
    
    # [SAFETY CHECK] Check raw occurrence in region before any heavy processing
    min_raw_obs <- 50
    
    if (nrow(coords_region) < min_raw_obs) {
      message(paste("SKIP: Not enough raw points for", sp_name, "in", region,
                    "(Found:", nrow(coords_region), "- Minimum required:", min_raw_obs, ")"))
      # Skips to next species in loop
      next
    }
    
    cat("Records in ", region, ":", nrow(coords_region), "\n")
    

    # -2.3- Remove occurrence points in the ocean --------------------------------------------------------------
    
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
    
    
    # -3- Environmental data and WorldClim bioclimatic variables -----------------------------------------------
    
    # -3.1- Define a study extent with a buffer ----------------------------------------------------------------
    
    region_ext <- ext(bounds[1], bounds[2], bounds[3], bounds[4])
    bioclim_crop <- crop(bioclim_global, region_ext)
    
    
    # -3.2- Extract raster values at occurrence points ---------------------------------------------------------
    
    # Ensure point CRS matches rasters
    species_pts <- project(species_land_vect, crs(bioclim_crop))
    
    # Extract climate values at points
    clim_vals <- terra::extract(bioclim_crop, species_pts)
    
    # Combine coordinates with climate values
    species_data <- bind_cols(
      as_tibble(crds(species_pts)) %>% rename(lon = x, lat = y),
      as_tibble(clim_vals)[,-1]
    ) %>% drop_na()
    
    # [SAFETY CHECK] Check processed occurrence points in region before generating current SDM
    min_obs <- 50
    
    if (nrow(species_data) < min_obs) {
      message(paste("SKIP: Not enough processed occurrence points for", sp_name,
                    "in", region, ". Only", nrow(species_data), "valid points found."))
      next
    }
    
    cat("Final number of dataset rows:", nrow(species_data), "\n")
    
    # Save processed data
    write.csv(species_data, here("data", "processed", paste0(sp_filename, ".csv")))
    
    
    # -4- Current Species Distribution Model Generation --------------------------------------------------------
    
    # This section includes building the GLM by generating pseudo-absence points (as we
    # don't have the data for actual absences) and splitting the presence/absence data
    # into training and testing data sets. The GLM fitted using the training data, and
    # evaluated via AUC scores using the testing data. Using this model, we can map a
    # prediction for habitat suitability.
    
    
    # -4.1- Generate background/pseudo-absence points ----------------------------------------------------------
    
    # Dynamic sample size for background points (either 1000 or number of presences if > 1000)
    bg_n <- max(1000, nrow(species_data))
    
    # Sample random points from the first layer of the cropped climate data
    set.seed(123)
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
    
    cat("Background points after NAs removed:", nrow(background_data), "\n")
    
    # Format presence data (Presence = 1)
    presence_data <- species_data %>% mutate(presence = 1)
    
    
    # -4.2- Split data into training (70%) and testing (30%) data
    
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
    
    
    # -4.3- Fit GLM --------------------------------------------------------------------------------------------
    
    # Identify variables for this species defined by the user
    current_vars <- predictor_list[[sp_name]]
    if(is.null(current_vars)) stop(paste("No variables stored for", sp_name))
    
    message(paste("Fitting GLM with: ", paste(current_vars, collapse = ", ")))
    
    # Construct formula dynamically
    formula_str <- paste("presence ~ ", paste(current_vars, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    # Fit model using training data
    sdm_model <- glm(model_formula, data = train_data, family = binomial)
    
    
    # -4.4- Evaluate model -------------------------------------------------------------------------------------
    
    # Evaluate using held-out test data
    eval_res <- evaluate(
      p = test_pres %>% dplyr::select(all_of(current_vars)),
      a = test_bg %>% dplyr::select(all_of(current_vars)),
      model = sdm_model
    )
    cat("AUC score:", eval_res@auc, "\n")
    
    
    # -4.5- Predict and Map ------------------------------------------------------------------------------------
    
    # Ensure raster has those layers
    if (!all(current_vars %in% names(bioclim_crop))) {
      stop("ERROR: One or more predictor layers not found in bioclim_crop: ", paste(current_vars, collapse = ", "))
    }
    
    # Predict suitability across the full study extent
    prediction <- terra::predict(bioclim_crop[[current_vars]], sdm_model, type = "response")
    
    # Store in list to contain map using the species name as a label
    current_sdm_maps[[sp_name]] <- prediction
    
    # Save stats to dataframe
    results_stats <- rbind(results_stats, data.frame(
      Species = sp_name,
      Region = region,
      N_Points = nrow(species_data),
      AUC = eval_res@auc
    ))
    
    # Plot
    plot(prediction, main = paste("SDM:", sp_name))
    points(species_land_vect, pch = 16, cex = 0.5, col = "black")
    
    message("Completed generation of current SDM for ", sp_name, " in ", region)
    
    cat("---------------------------------------------------------------------------------------------------\n")

  }
  
  message("All SDMs generated successfully. Returning map list...")
  return(list(maps = current_sdm_maps, stats = results_stats))
  
}


# Run the function
output <- run_current_sdm(sp1, sp2, "Europe", sp_predictors)

# Get the maps for Task 2
my_maps <- output$maps

# Get the stats table for the current SDMs
print(output$stats)






































































