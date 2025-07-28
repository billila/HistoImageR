#' Import HoVerNet output (JSON) as SpatialExperiment
#' 
#' @importFrom jsonlite fromJSON
#' @importFrom readr read_file
#' 
#' @param fnames A character vector containing file paths to image feature files.
#' This function expects that file names include the TCGA patient ID in a `TCGA-\\d{2}-\\d{4}` pattern.   
#' @param subset A numeric value >= 1. `1/subset` amount of cells will be collected.
#' If you want to import all the cells, set this as `1`.
#' 
#' @return A SpatialExperiment object
#' 
#'  
JsonToSe <- function(fnames, subset = subset) {
  
  text_content <- readr::read_file(fnames)
  json_data <- jsonlite::fromJSON(text_content)
  cell_list <- json_data[[2]] # a list of nucleus - `nuc` slot from HoVerNet output
  
  randomSubset <- sample(seq_along(cell_list), round(length(cell_list)/subset, 0)) # Use only a subset of cells for demo (default: 1/10)
  x <- cell_list[sort(randomSubset)]
  
  ## Create SpatialExperiment object
  se <- SpatialExperiment(
    spatialCoords = do.call(rbind, lapply(x, function(x) x$centroid)), # 1. coordinates for centroid
    colData = DataFrame(
      cell_type = sapply(x, function(x) x$type), # 2. cell type
      type_prob = sapply(x, function(x) x$type_prob) # 3. cell type probability
    )
  )
  
  ## Reformat cell labels/numbers
  colnames(se) <- 
  
  ## Store segmentation polygons in metadata
  metadata(se)$bboxes <- lapply(x, function(x) x$bbox) # 4. bbox
  metadata(se)$polygons <- lapply(x, function(x) x$contour) # 5. contour
  
  ## Update column names
  colnames(spatialCoords(se)) <- c("x", "y") # spatialCoords slot
  
  ## Convert cell_type as factors
  colData(se)$cell_type <- as.factor(colData(se)$cell_type) 
  
  return(se)
}