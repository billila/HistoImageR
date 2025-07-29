#' Import Squidpy output (h5ad) as SpatialExperiment
#'
#' @importFrom anndata read_h5ad
#'
#' @param fname A character(1). Containing a file path to an image feature file.
#' This function expects that file names include the TCGA patient ID in a `TCGA-\\d{2}-\\d{4}` pattern.   
#' @param subsample_rate A numeric value >= 1. `1/subsample_rate` amount of cells will be collected.
#' If you want to import all the cells, set this as `1`.
#' @param seed Integer or NULL. Random seed for reproducible results. 
#' If NULL (default), no seed is set and results will be random.
#' 
#' 
#' 
h5adToSe <- function(fname, subsample_rate = subsample_rate, seed = seed) {
    
    ann <- anndata::read_h5ad(fname)
    
    ## Create SpatialExperiment object
    spe <- SpatialExperiment(
        spatialCoords = ann$obsm$spatial, # 1. coordinates for centroid
        colData = ann$obs # 2. cell type
    )
    
    # ## Store segmentation polygons in metadata
    # metadata(spe)$bboxes <- lapply(ann, function(x) x$bbox) # 4. bbox
    # metadata(spe)$polygons <- lapply(ann, function(x) x$contour) # 5. contour
    
    ## Update column names
    colnames(spatialCoords(spe)) <- c("x", "y") # spatialCoords slot
    colnames(colData(spe))[1] <- "cell_type"
    
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