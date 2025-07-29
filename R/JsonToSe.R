#' Import HoVerNet output (JSON) as SpatialExperiment
#' 
#' @importFrom jsonlite fromJSON
#' @importFrom readr read_file
#' 
#' @param fname A character(1). Containing a file path to an image feature file.
#' This function expects that file names include the TCGA patient ID in a `TCGA-\\d{2}-\\d{4}` pattern.   
#' @param subsample_rate A numeric value >= 1. `1/subsample_rate` amount of cells will be collected.
#' If you want to import all the cells, set this as `1`.
#' @param seed Integer or NULL. Random seed for reproducible results. 
#' If NULL (default), no seed is set and results will be random.
#' 
#' @return A SpatialExperiment object
#' 
#' @export
jsonToSe <- function(fname, subsample_rate = subsample_rate, seed = seed) {
    
    ## Load JSON file as a list
    text_content <- readr::read_file(fname)
    json_data <- jsonlite::fromJSON(text_content)
    cell_list <- json_data[["nuc"]] # a list of nucleus
    names(cell_list) <- seq_along(cell_list) # reformat cell labels/numbers for consistency
  
    ## Create SpatialExperiment object
    spe <- SpatialExperiment(
        spatialCoords = do.call(rbind, lapply(cell_list, function(x) x$centroid)), # 1. coordinates for centroid
        colData = DataFrame(
            cell_type = sapply(cell_list, function(x) x$type), # 2. cell type
            type_prob = sapply(cell_list, function(x) x$type_prob) # 3. cell type probability
        )
    )
  
    ## Store segmentation polygons in metadata
    metadata(spe)$bboxes <- lapply(cell_list, function(x) x$bbox) # 4. bbox
    metadata(spe)$polygons <- lapply(cell_list, function(x) x$contour) # 5. contour
  
    ## Update column names
    colnames(spatialCoords(spe)) <- c("x", "y") # spatialCoords slot
  
    ## Convert cell_type as factors
    colData(spe)$cell_type <- as.factor(colData(spe)$cell_type) 
    
    ## Subset for quicker processing
    if (!is.null(seed)) {set.seed(seed)} # seed for random subsetting
    
    ####### anndata object indexing needs to be matched with json <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    colnames(spe) <- seq_len(ncol(spe))
    sample_size <- floor(ncol(spe) / subsample_rate)
    randomSubset <- sample(ncol(spe), sample_size, replace = FALSE) # Use only a subset of cells for demo (default: 1/10)
    spe <- spe[,sort(randomSubset)]
  
    return(spe)
}