#' Import HoVerNet output (JSON) as SpatialFeatureExperiment
#' 
#' @importFrom jsonlite fromJSON
#' @importFrom readr read_file
#' @importFrom dplyr rename
#' 
#' @param fname A character(1). Containing a file path to an image feature file.
#' This function expects that file names include the TCGA patient ID in a `TCGA-\\d{2}-\\d{4}` pattern.   
#' @param subsample_rate A numeric value >= 1. `1/subsample_rate` amount of cells will be collected.
#' If you want to import all the cells, set this as `1`.
#' @param seed Integer or NULL. Random seed for reproducible results. 
#' If NULL (default), no seed is set and results will be random.
#' 
#' @return A SpatialFeatureExperiment object
#' 
#'  
jsonToSfe <- function(fname, subsample_rate = subsample_rate, seed = seed) {
    
    text_content <- readr::read_file(fname)
    json_data <- jsonlite::fromJSON(text_content)
    cell_list <- json_data[["nuc"]] # a list of nucleus
    names(cell_list) <- seq_along(cell_list) # reformat cell labels/numbers for consistency

    ## Adding centroids (as sf POINT objects)
    cell_coords <- t(sapply(cell_list, function(x) x$centroid)) %>%
        as.data.frame(.) %>%
        dplyr::rename(., x = V1, y = V2)

    ## Adding bounding boxes (as sf POLYGON objects)
    bbox_polygons <- lapply(cell_list, function(x) {
        bbox_coords <- x$bbox
        xmin <- bbox_coords[1,1]
        ymin <- bbox_coords[1,2]
        xmax <- bbox_coords[2,1]
        ymax <- bbox_coords[2,2]
        st_polygon(list(matrix(c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin), ncol = 2, byrow = TRUE)))
    })
    bbox_polygons <- st_sfc(bbox_polygons)
    
    ## Cell segmentation polygons (contours)
    contour_polygons <- lapply(cell_list, function(x) {
        contour_polygon_vector <- x$contour %>% t(.) %>% as.vector(.)
        first_point <- head(contour_polygon_vector, 2) # extract the first point
        closed_contour_polygon_vector <- c(contour_polygon_vector, first_point) # append the first point to close the polygon
        st_polygon(list(matrix(closed_contour_polygon_vector, ncol = 2, byrow = TRUE)))
    })
    contour_polygons <- st_sfc(contour_polygons)
    
    ## Empty matrix as a placeholder
    empty_matrix <- matrix(0, nrow = 1, ncol = length(cell_list))
    colnames(empty_matrix) <- seq_along(cell_list)
    
    ## Create SpatialExperiment object
    sfe <- SpatialFeatureExperiment(
        assays = list(assay = empty_matrix), 
        spatialCoords = do.call(rbind, lapply(cell_list, function(x) x$centroid)), # 1. coordinates for centroid
        colData = DataFrame(
            cell_type = sapply(cell_list, function(x) x$type), # 2. cell type
            type_prob = sapply(cell_list, function(x) x$type_prob) # 3. cell type probability
        )
    )
    
    ## Add the spatial features as colGeometries
    sfe$colGeometries$cell_centroids <- st_as_sf(cell_coords, coords = c("x", "y")) # <<<<<<<<<<<<<<<<<<<< this takes really long
    sfe$colGeometries$cell_bboxes <- st_as_sf(bbox_polygons, crs = st_crs(sfe$colGeometries$cell_centroids)) # bbox <<<<<<<<<<<<<<<<<<<<<<< BUG here
    sfe$colGeometries$cell_contours <- st_as_sf(contour_polygons, crs = st_crs(sfe$colGeometries$cell_centroids)) # contour
    
    ## Adding a tissue boundary polygon <<<<<<<<<<<<<< we don't have this yet
    # sfe$annotGeometries$tissue_boundary <- tissue_polygon
  
    ## Convert cell_type as factors
    colData(spe)$cell_type <- as.factor(colData(spe)$cell_type) 
  
    return(spe)
}