# Nuclei Feature Extraction from H&E Images using R/Bioconductor
# Updated for nucleus-wise input structure
# 
# Author: Feature Extraction Pipeline
# Date: 2025
# 
# Description: Comprehensive feature extraction from H&E stained nuclei data
# Input: List of n nuclei, each containing bbox, centroid, contour, type_prob, type
# Output: Extracted morphological, spatial, and classification features

# Required packages
suppressPackageStartupMessages({
    library(EBImage)      # Image processing and morphological measurements
    library(spatstat)     # Spatial point pattern analysis
    library(sf)           # Simple features for geometric computations
    library(dplyr)        # Data manipulation
    library(moments)      # Statistical moments for shape analysis
    library(pracma)       # Practical numerical math functions
})

#' Extract comprehensive features from nuclei data in H&E images
#'
#' This function processes a list of nuclei, each containing bounding box, centroid,
#' contour, type probability, and type information, and extracts morphological,
#' spatial, and classification features useful for digital pathology analysis.
#'
#' @param nuclei_list A list of length n, where each element represents one nucleus
#'                    containing the following named elements:
#'                    - bbox: numeric vector c(xmin, ymin, xmax, ymax)
#'                    - centroid: numeric vector c(x, y)
#'                    - contour: matrix with columns x, y representing boundary points
#'                    - type_prob: numeric value (0-1) indicating classification confidence
#'                    - type: character string indicating nucleus type
#' @param image_dims Numeric vector c(width, height) specifying image dimensions
#'                   for spatial analysis. Default: c(1024, 1024)
#' @param verbose Logical indicating whether to print progress messages. Default: TRUE
#'
#' @return A list containing extracted features organized by category:
#'         - individual: data.frame with per-nucleus features
#'         - spatial: list with population-level spatial features
#'         - population: list with population-level statistics
#'
#' @examples
#' # Create example nucleus
#' nucleus1 <- list(
#'   bbox = c(100, 150, 130, 180),
#'   centroid = c(115, 165),
#'   contour = matrix(c(seq(110,120,length=20), seq(160,170,length=20)), ncol=2),
#'   type_prob = 0.85,
#'   type = "epithelial"
#' )
#' nuclei_data <- list(nucleus1)  # Can contain multiple nuclei
#' features <- extract_nuclei_features(nuclei_data)
#'
extract_nuclei_features <- function(nuclei_list, 
                                    image_dims = c(1024, 1024), 
                                    verbose = TRUE) {
    
    if(verbose) cat("Starting feature extraction for", length(nuclei_list), "nuclei...\n")
    
    # Validate input
    if(!is.list(nuclei_list) || length(nuclei_list) == 0) {
        stop("nuclei_list must be a non-empty list")
    }
    
    # Check required fields in first nucleus
    required_fields <- c("bbox", "centroid", "contour", "type_prob", "type")
    if(!all(required_fields %in% names(nuclei_list[[1]]))) {
        missing_fields <- setdiff(required_fields, names(nuclei_list[[1]]))
        stop("Missing required fields in nuclei: ", paste(missing_fields, collapse = ", "))
    }
    
    n_nuclei <- length(nuclei_list)
    
    # Initialize results
    features <- list(
        individual = data.frame(),
        spatial = list(),
        population = list()
    )
    
    # ====== INDIVIDUAL NUCLEUS FEATURES ======
    if(verbose) cat("Extracting individual nucleus features...\n")
    
    individual_features <- extract_individual_features(nuclei_list, verbose)
    features$individual <- individual_features
    
    # ====== SPATIAL FEATURES ======
    if(verbose) cat("Computing spatial pattern features...\n")
    
    spatial_features <- extract_spatial_features(nuclei_list, image_dims, verbose)
    features$spatial <- spatial_features
    
    # ====== POPULATION FEATURES ======
    if(verbose) cat("Computing population-level features...\n")
    
    population_features <- extract_population_features(nuclei_list, verbose)
    features$population <- population_features
    
    if(verbose) cat("Feature extraction completed successfully!\n")
    
    return(features)
}

#' Extract features for individual nuclei
#'
#' Computes morphological and classification features for each nucleus separately.
#'
#' @param nuclei_list List of nuclei with required fields
#' @param verbose Logical for progress messages
#'
#' @return data.frame with one row per nucleus containing individual features
#'
#' @details Extracted features include:
#'          - Bounding box metrics (area, aspect ratio, etc.)
#'          - Shape descriptors from contours (circularity, solidity, etc.)
#'          - Classification confidence metrics
extract_individual_features <- function(nuclei_list, verbose = FALSE) {
    
    n_nuclei <- length(nuclei_list)
    
    # Initialize feature matrix
    individual_df <- data.frame(
        nucleus_id = 1:n_nuclei,
        # Bounding box features
        bbox_width = numeric(n_nuclei),
        bbox_height = numeric(n_nuclei),
        bbox_area = numeric(n_nuclei),
        bbox_aspect_ratio = numeric(n_nuclei),
        bbox_perimeter = numeric(n_nuclei),
        
        # Contour-based shape features
        contour_area = numeric(n_nuclei),
        contour_perimeter = numeric(n_nuclei),
        circularity = numeric(n_nuclei),
        solidity = numeric(n_nuclei),
        convexity = numeric(n_nuclei),
        eccentricity = numeric(n_nuclei),
        major_axis = numeric(n_nuclei),
        minor_axis = numeric(n_nuclei),
        elongation = numeric(n_nuclei),
        compactness = numeric(n_nuclei),
        
        # Classification features
        type_prob = numeric(n_nuclei),
        type = character(n_nuclei),
        
        stringsAsFactors = FALSE
    )
    
    for(i in 1:n_nuclei) {
        nucleus <- nuclei_list[[i]]
        
        # Extract bounding box features
        bbox_features <- compute_bbox_features(nucleus$bbox)
        individual_df[i, names(bbox_features)] <- bbox_features
        
        # Extract contour features
        contour_features <- compute_contour_features(nucleus$contour)
        individual_df[i, names(contour_features)] <- contour_features
        
        # Store classification info
        individual_df$type_prob[i] <- nucleus$type_prob
        individual_df$type[i] <- nucleus$type
        
        if(verbose && i %% 50 == 0) {
            cat("  Processed", i, "of", n_nuclei, "nuclei\n")
        }
    }
    
    return(individual_df)
}

#' Compute bounding box-based features
#'
#' @param bbox Numeric vector c(xmin, ymin, xmax, ymax)
#' @return Named list of bounding box features
compute_bbox_features <- function(bbox) {
    
    width <- bbox[3] - bbox[1]
    height <- bbox[4] - bbox[2]
    area <- width * height
    perimeter <- 2 * (width + height)
    aspect_ratio <- width / height
    
    return(list(
        bbox_width = width,
        bbox_height = height,
        bbox_area = area,
        bbox_aspect_ratio = aspect_ratio,
        bbox_perimeter = perimeter
    ))
}

#' Compute contour-based morphological features
#'
#' @param contour Matrix with x,y coordinates of nucleus boundary
#' @return Named list of shape features
#'
#' @details Computes comprehensive shape descriptors including:
#'          - Basic metrics (area, perimeter)
#'          - Shape indices (circularity, solidity, convexity)
#'          - Ellipse fitting parameters (eccentricity, axes)
compute_contour_features <- function(contour) {
    
    # Handle edge cases
    if(nrow(contour) < 3) {
        return(list(
            contour_area = 0, contour_perimeter = 0, circularity = 0,
            solidity = 0, convexity = 0, eccentricity = 0,
            major_axis = 0, minor_axis = 0, elongation = 0, compactness = 0
        ))
    }
    
    # Basic geometric measurements
    area <- tryCatch({
        abs(polyarea(contour[,1], contour[,2]))
    }, error = function(e) {
        # If polyarea fails, estimate area using bounding box
        x_range <- diff(range(contour[,1], na.rm = TRUE))
        y_range <- diff(range(contour[,2], na.rm = TRUE))
        x_range * y_range * 0.785  # Approximate ellipse area
    })
    
    # Handle NA or invalid area values
    if(is.na(area) || is.null(area) || !is.finite(area)) {
        area <- 1e-10
    }
    if(area <= 0) area <- 1e-10
    
    # Perimeter calculation
    dx <- diff(c(contour[,1], contour[1,1]))
    dy <- diff(c(contour[,2], contour[1,2]))
    perimeter <- sum(sqrt(dx^2 + dy^2), na.rm = TRUE)
    
    # Handle invalid perimeter values
    if(is.na(perimeter) || is.null(perimeter) || !is.finite(perimeter) || perimeter <= 0) {
        perimeter <- 1e-10
    }
    
    # Circularity: 4π*Area/Perimeter²
    circularity <- (4 * pi * area) / (perimeter^2)
    
    # Convex hull analysis
    solidity <- 1.0
    convexity <- 1.0
    
    tryCatch({
        # Check for valid contour data
        if(nrow(contour) >= 3 && !any(is.na(contour))) {
            hull_indices <- chull(contour)
            if(length(hull_indices) >= 3) {
                hull_contour <- contour[hull_indices, ]
                hull_area <- abs(polyarea(hull_contour[,1], hull_contour[,2]))
                hull_perimeter <- sum(sqrt(diff(c(hull_contour[,1], hull_contour[1,1]))^2 + 
                                               diff(c(hull_contour[,2], hull_contour[1,2]))^2), na.rm = TRUE)
                
                # Check for valid hull calculations
                if(!is.na(hull_area) && is.finite(hull_area) && hull_area > 0) {
                    solidity <- area / hull_area
                }
                
                if(!is.na(hull_perimeter) && is.finite(hull_perimeter) && hull_perimeter > 0) {
                    convexity <- hull_perimeter / perimeter
                }
            }
        }
    }, error = function(e) {
        # Keep default values if convex hull calculation fails
    })
    
    # Ensure valid values
    if(is.na(solidity) || !is.finite(solidity) || solidity <= 0) solidity <- 1.0
    if(is.na(convexity) || !is.finite(convexity) || convexity <= 0) convexity <- 1.0
    
    # Ellipse fitting for eccentricity and axes
    ellipse_params <- fit_ellipse_to_contour(contour)
    
    # Elongation (aspect ratio)
    elongation <- 1.0
    if(!is.na(ellipse_params$minor_axis) && is.finite(ellipse_params$minor_axis) && ellipse_params$minor_axis > 0) {
        elongation <- ellipse_params$major_axis / ellipse_params$minor_axis
    }
    
    # Compactness
    compactness <- (perimeter^2) / (4 * pi * area)
    if(is.na(compactness) || !is.finite(compactness)) compactness <- 1.0
    
    return(list(
        contour_area = area,
        contour_perimeter = perimeter,
        circularity = circularity,
        solidity = solidity,
        convexity = convexity,
        eccentricity = ellipse_params$eccentricity,
        major_axis = ellipse_params$major_axis,
        minor_axis = ellipse_params$minor_axis,
        elongation = elongation,
        compactness = compactness
    ))
}

#' Fit ellipse to contour and extract parameters
#'
#' @param contour Matrix of x,y boundary coordinates
#' @return List with ellipse parameters
fit_ellipse_to_contour <- function(contour) {
    
    tryCatch({
        # Calculate centroid
        cx <- mean(contour[,1], na.rm = TRUE)
        cy <- mean(contour[,2], na.rm = TRUE)
        
        # Check for valid centroid
        if(is.na(cx) || is.na(cy)) {
            return(list(major_axis = 1.0, minor_axis = 1.0, eccentricity = 0.0))
        }
        
        # Calculate second moments
        x_centered <- contour[,1] - cx
        y_centered <- contour[,2] - cy
        
        # Remove any NA values
        valid_idx <- !is.na(x_centered) & !is.na(y_centered)
        if(sum(valid_idx) < 3) {
            return(list(major_axis = 1.0, minor_axis = 1.0, eccentricity = 0.0))
        }
        
        x_centered <- x_centered[valid_idx]
        y_centered <- y_centered[valid_idx]
        
        mu20 <- mean(x_centered^2)
        mu02 <- mean(y_centered^2)
        mu11 <- mean(x_centered * y_centered)
        
        # Check for valid moments
        if(any(is.na(c(mu20, mu02, mu11))) || any(!is.finite(c(mu20, mu02, mu11)))) {
            return(list(major_axis = 1.0, minor_axis = 1.0, eccentricity = 0.0))
        }
        
        # Eigenvalue analysis for ellipse parameters
        term1 <- (mu20 + mu02) / 2
        term2 <- sqrt(((mu20 - mu02) / 2)^2 + mu11^2)
        
        lambda1 <- term1 + term2
        lambda2 <- term1 - term2
        
        # Prevent negative or zero eigenvalues
        lambda1 <- max(lambda1, 1e-10)
        lambda2 <- max(lambda2, 1e-10)
        
        major_axis <- 2 * sqrt(lambda1)
        minor_axis <- 2 * sqrt(lambda2)
        
        # Eccentricity
        eccentricity <- sqrt(1 - (lambda2 / lambda1))
        
        # Ensure all values are finite
        if(any(!is.finite(c(major_axis, minor_axis, eccentricity)))) {
            return(list(major_axis = 1.0, minor_axis = 1.0, eccentricity = 0.0))
        }
        
        return(list(
            major_axis = major_axis,
            minor_axis = minor_axis,
            eccentricity = eccentricity
        ))
        
    }, error = function(e) {
        return(list(
            major_axis = 1.0,
            minor_axis = 1.0,
            eccentricity = 0.0
        ))
    })
}

#' Extract spatial pattern features from nucleus centroids
#'
#' @param nuclei_list List of nuclei
#' @param image_dims Image dimensions c(width, height)
#' @param verbose Progress messages flag
#'
#' @return List of spatial analysis results
#'
#' @details Performs spatial point pattern analysis including:
#'          - Nuclear density calculations
#'          - Nearest neighbor analysis
#'          - Clustering analysis using Ripley's K-function
#'          - Spatial regularity metrics
extract_spatial_features <- function(nuclei_list, image_dims, verbose = FALSE) {
    
    # Extract centroids
    centroids <- t(sapply(nuclei_list, function(n) n$centroid))
    
    # Create spatial point pattern
    x_coords <- centroids[,1]
    y_coords <- centroids[,2]
    
    # Define observation window
    W <- owin(c(0, image_dims[1]), c(0, image_dims[2]))
    
    # Create point pattern object
    pp <- ppp(x_coords, y_coords, window = W)
    
    spatial_results <- list()
    
    # Basic density metrics
    spatial_results$nuclear_density <- pp$n / area(W)
    spatial_results$n_nuclei <- pp$n
    
    # Nearest neighbor analysis
    if(pp$n > 1) {
        nn_distances <- nndist(pp)
        spatial_results$mean_nn_distance <- mean(nn_distances)
        spatial_results$median_nn_distance <- median(nn_distances)
        spatial_results$sd_nn_distance <- sd(nn_distances)
        spatial_results$min_nn_distance <- min(nn_distances)
        spatial_results$max_nn_distance <- max(nn_distances)
        
        # Regularity index (Clark-Evans)
        expected_nn <- 1 / (2 * sqrt(spatial_results$nuclear_density))
        spatial_results$clark_evans_index <- spatial_results$mean_nn_distance / expected_nn
        
        # Spatial clustering analysis
        if(pp$n >= 10) {  # Need sufficient points for K-function
            tryCatch({
                K_est <- Kest(pp, correction = "isotropic")
                # Maximum deviation from complete spatial randomness
                spatial_results$ripley_k_max_deviation <- max(K_est$iso - K_est$theo, na.rm = TRUE)
                spatial_results$ripley_k_clustering <- spatial_results$ripley_k_max_deviation > 0
                
            }, error = function(e) {
                if(verbose) cat("  Warning: Could not compute Ripley's K function\n")
                spatial_results$ripley_k_max_deviation <<- NA
                spatial_results$ripley_k_clustering <<- NA
            })
        }
    } else {
        # Single point case
        spatial_results$mean_nn_distance <- NA
        spatial_results$median_nn_distance <- NA
        spatial_results$sd_nn_distance <- NA
        spatial_results$clark_evans_index <- NA
    }
    
    return(spatial_results)
}

#' Extract population-level features from all nuclei
#'
#' @param nuclei_list List of nuclei
#' @param verbose Progress messages flag
#'
#' @return List of population statistics
extract_population_features <- function(nuclei_list, verbose = FALSE) {
    
    # Extract type information
    types <- sapply(nuclei_list, function(n) n$type)
    type_probs <- sapply(nuclei_list, function(n) n$type_prob)
    
    population_results <- list()
    
    # Classification confidence statistics
    population_results$mean_confidence <- mean(type_probs)
    population_results$median_confidence <- median(type_probs)
    population_results$min_confidence <- min(type_probs)
    population_results$max_confidence <- max(type_probs)
    population_results$confidence_std <- sd(type_probs)
    
    # Type distribution analysis
    type_counts <- table(types)
    population_results$n_cell_types <- length(type_counts)
    population_results$dominant_type <- names(type_counts)[which.max(type_counts)]
    population_results$dominant_type_fraction <- max(type_counts) / sum(type_counts)
    
    # Shannon diversity index for cell types
    type_proportions <- type_counts / sum(type_counts)
    population_results$shannon_diversity <- -sum(type_proportions * log(type_proportions + 1e-10))
    
    # Simpson diversity index
    population_results$simpson_diversity <- 1 - sum(type_proportions^2)
    
    # Type-specific confidence analysis
    type_confidence_stats <- aggregate(type_probs, by = list(Type = types), 
                                       FUN = function(x) c(mean = mean(x), std = sd(x)))
    population_results$type_confidence_stats <- type_confidence_stats
    
    return(population_results)
}

#' Create example nuclei data for testing
#'
#' @param n_nuclei Number of nuclei to generate
#' @param image_dims Image dimensions for realistic positioning
#' @param seed Random seed for reproducibility
#'
#' @return List of example nuclei with required structure
create_example_nuclei <- function(n_nuclei = 100, image_dims = c(1024, 1024), seed = 123) {
    
    set.seed(seed)
    
    nuclei_list <- list()
    
    cell_types <- c("epithelial", "stromal", "immune", "endothelial")
    type_probs <- c(0.4, 0.3, 0.2, 0.1)
    
    for(i in 1:n_nuclei) {
        # Generate random position
        center_x <- runif(1, 50, image_dims[1] - 50)
        center_y <- runif(1, 50, image_dims[2] - 50)
        
        # Generate size parameters
        semi_major <- runif(1, 8, 20)
        semi_minor <- runif(1, 6, semi_major)
        
        # Create bounding box
        bbox <- c(center_x - semi_major, center_y - semi_minor,
                  center_x + semi_major, center_y + semi_minor)
        
        # Create contour (elliptical)
        theta <- seq(0, 2*pi, length.out = 50)
        contour_x <- center_x + semi_major * cos(theta) + rnorm(50, 0, 0.5)
        contour_y <- center_y + semi_minor * sin(theta) + rnorm(50, 0, 0.5)
        contour <- cbind(contour_x, contour_y)
        
        # Assign cell type and confidence
        cell_type <- sample(cell_types, 1, prob = type_probs)
        confidence <- runif(1, 0.6, 0.99)
        
        # Create nucleus object
        nuclei_list[[i]] <- list(
            bbox = bbox,
            centroid = c(center_x, center_y),
            contour = contour,
            type_prob = confidence,
            type = cell_type
        )
    }
    
    return(nuclei_list)
}







# # ====== DEMONSTRATION AND TESTING ======
# 
# cat("=== NUCLEI FEATURE EXTRACTION PIPELINE ===\n")
# cat("Creating example data...\n")
# 
# # Generate example data
# example_nuclei <- create_example_nuclei(n_nuclei = 150, seed = 42)
# 
# cat("Example nucleus structure:\n")
# str(example_nuclei[[1]], max.level = 2)
# 
# # Extract features
# cat("\n=== RUNNING FEATURE EXTRACTION ===\n")
# example_nuclei <- json_data[[2]]
# extracted_features <- extract_nuclei_features(example_nuclei, 
#                                               image_dims = c(1024, 1024),
#                                               verbose = TRUE)
# 
# # Display results summary
# cat("\n=== RESULTS SUMMARY ===\n")
# 
# cat("\nIndividual Features (first 5 nuclei):\n")
# print(head(extracted_features$individual[, 1:8], 5))
# 
# cat("\nSpatial Features:\n")
# # Handle potential NA values in spatial features
# spatial_values <- sapply(extracted_features$spatial, function(x) {
#     if(is.numeric(x) && !is.na(x)) {
#         return(as.character(round(x, 4)))
#     } else {
#         return(as.character(x))
#     }
# })
# 
# spatial_df <- data.frame(
#     Feature = names(extracted_features$spatial),
#     Value = spatial_values,
#     stringsAsFactors = FALSE
# )
# print(spatial_df)
# 
# cat("\nPopulation Features:\n")
# pop_simple <- extracted_features$population[!names(extracted_features$population) %in% "type_confidence_stats"]
# 
# # Handle mixed data types (numeric and character)
# pop_values <- sapply(pop_simple, function(x) {
#     if(is.numeric(x)) {
#         return(as.character(round(x, 4)))
#     } else {
#         return(as.character(x))
#     }
# })
# 
# pop_df <- data.frame(
#     Feature = names(pop_simple),
#     Value = pop_values,
#     stringsAsFactors = FALSE
# )
# print(pop_df)
# 
# cat("\n=== FEATURE EXTRACTION PIPELINE READY FOR USE ===\n")
# cat("Total features extracted per nucleus:", ncol(extracted_features$individual) - 1, "\n")
# cat("Use extract_nuclei_features(your_nuclei_list) with your data!\n")