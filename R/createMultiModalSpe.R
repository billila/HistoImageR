#' Merging spatially-resolved image features and sample-level data
#' 
#' This package combine multi-omics data for TCGA samples and extracted 
#' TCGA image features into a `SpatialFeatureExperiment` object.
#' 
#' FYI, DNA methylation data comes closest to having the most comprehensive 
#' coverage across TCGA samples with diagnostic images.
#'  
#' @param x An `ExperimentList` objects including `SpatialExperiment` objects. 
#' @param meta A data frame with participants' clinical metadata
#' 
#' @examples
#' # ovmae <- curatedTCGAData::curatedTCGAData("OV", "Methylation*", version = "2.0.1", dry.run = FALSE)
#' # meta <- colData(ovmae)
#' 
createMultiModalSpe <- function(x, meta) {
    
    ## Pre-allocate a list to store updated SpatialExperiment objects
    updated_x <- vector("list", length(x))
    names(updated_x) <- names(x)
    
    ## Add patientID to cell_ids
    for (i in seq_along(x)) {
        
        patient_id <- names(x)[i]
        spe <- x[[i]] 
        
        ## Organize image metadata        
        df <- colData(spe)
        df$cell_number <- rownames(df)
        df$patientID <- patient_id
        df$cell_id <- paste(patient_id, colnames(spe), sep = "_")
        
        ## Merge with clinical metadata
        mergedMeta <- merge(df, meta, by = "patientID", all.x = TRUE)
        
        ## Update colData
        colData(spe) <- mergedMeta
        colnames(spe) <- mergedMeta$cell_id
        
        ## Store multi-omics assay data into metadata slot <<<<<<<<<<
        
        ## Store updated object
        updated_x[[i]] <- spe
    }
    
    ## Combine all SpatialExperiment objects
    res <- do.call(cbind, updated_x)
    
    return(res)
}