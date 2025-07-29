#' Analyze the characteristics of embedding clusters
#' 
#' 
#' 
#' 
analyze_embedding_separation <- function(embedding_matrix,
                                         metadata_df,
                                         metadata_cols = NULL) {
    # Check inputs
    if (nrow(embedding_matrix) != nrow(metadata_df)) {
        stop("Number of samples in embedding matrix and metadata must match")
    }
    
    # If no specific columns specified, use all non-numeric columns
    if (is.null(metadata_cols)) {
        metadata_cols <- names(metadata_df)[sapply(metadata_df, function(x) !is.numeric(x))]
    }
    
    results <- list()
    
    # PCA
    pca_result <- prcomp(embedding_matrix, center = TRUE, scale. = TRUE)
    pca_df <- data.frame(PC1 = pca_result$x[, 1], 
                         PC2 = pca_result$x[, 2],
                         metadata_df)
    
    plots <- list()
    
    for (col in metadata_cols) {
        if (length(unique(metadata_df[[col]])) > 1) {
            p <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = col)) +
                geom_point(alpha = 0.7, size = 2) +
                theme_minimal() +
                labs(title = paste("PCA visualization colored by", col),
                     x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
                     y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
                theme(legend.position = "bottom")
            
            plots[[paste0("pca_", col)]] <- p
        }
    }
    
    results$plots <- plots
    results$pca_result <- pca_result
    
    return(results)
}