#' Merging spatially-resolved image features and sample-level data
#' 
#' This package combine multi-omics data for TCGA samples and extracted 
#' TCGA image features into a `SpatialFeatureExperiment` object.
#' 
#' FYI, DNA methylation data comes closest to having the most comprehensive 
#' coverage across TCGA samples with diagnostic images.
#'  
#' 
#' @param x An `ExperimentList` objects including `SpatialExperiment` objects. 
#' @param meta A data frame with participants' clinical metadata
#' 
createMultiModalSpe <- function(x, meta) {
  
  ## Add patientID to cell_ids
  for (i in seq_along(x)) {
    
    ## Organize image metadata        
    df <- colData(x[[i]])
    df$cell_number <- rownames(df)
    cellID <- paste(names(x)[i], colnames(x[[i]]), sep = "_")
    df$cell_id <- cellID
    df$patientID <- names(x)[i]
    
    ## Merge with clinical metadata
    mergedMeta <- merge(df, meta, by = "patientID", all.x = TRUE)
    
    ## Update colData
    colData(x[[i]]) <- mergedMeta
    colnames(x[[i]]) <- cellID
    
    ## Store multi-omics assay data into metadata slot
    
  }
  
  ## Combine all SpatialExperiment objects
  res <- do.call(cbind, x)
  
  return(res)
}