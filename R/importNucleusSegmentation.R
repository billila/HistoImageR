.importH5ad <- function(fname) {
  ann <- anndata::read_h5ad(h5ad_file)
  return(ann)
}

.importJson <- function(fname, subset = subset) {
  text_content <- readr::read_file(fnames)
  json_data <- jsonlite::fromJSON(text_content)
  cells <- json_data[[2]] # a list of nucleus
  
  randomSubset <- sample(seq_along(cells), round(length(cells)/subset, 0)) # Use only a subset of cells for demo (default: 1/10)
  spe <- cellListToSE(cells[sort(randomSubset)])
  
  return(spe)
}

#' Import nucleu features for TCGA samples as SpatialExperiment
#' 
#' This function imports HoVerNet outputs (JSON) or Squidpy outputs (h5ad) as R list.
#' 
#' @import jsonlite
#' @import SpatialExperiment
#' 
#' @param x A character vector containing file paths to image feature files.
#' @param subset A numeric value >= 1. `1/subset` amount of cells will be collected.
#' If you want to import all the cells, set this as `1`.
#' 
#' @return An ExperimentList object containing SpatialExperiment objects (one per slide/sample)
#' 
importNucleusSegmentation <- function(fnames, subset = 10) {
  
  patientIDs <- stringr::str_extract(fnames, "TCGA-\\d{2}-\\d{4}")
  spelist <- vector(mode = "list", length = length(fnames))
  names(spelist) <- patientIDs
  
  for (i in seq_along(fnames)) {
    
    dat <- fnames[i]
    
    if (grepl(".json$", dat)) {
      spe <- JsonToSe(dat)
    } else if (grepl(".h5ad$", dat)) {
      spe <- H5adToSe(dat)
    }

    message(paste(i, "out of", length(fnames), "image has imported"))
    
    # ## Create SotredSpatialImage object (from local file)
    # image_file <- file.path(pathToFiles, paste0(fnames[i], ".png"))
    # spe <- addImg(spe, 
    #               imageSource = image_file,
    #               sample_id = "sample01",
    #               image_id = "h&e",
    #               scaleFactor = 1)
    
    ## Combine all results as ExperimentList
    spelist[i] <- spe
  }
  
  spes <- ExperimentList(spelist)
  
}