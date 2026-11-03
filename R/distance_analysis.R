#' Calculate Distance Profile
#'
#' Calculates genomic separation metrics at varying distances from boundary
#'
#' @param seurat_obj Seurat object with spatial data
#' @param distance_thresholds Vector of distances to test (in microns)
#' @param metric Which metrics to calculate ("all", "transcriptomic", "composition", "de_strength")
#' @param group_by Metadata column for cell type composition
#' @param inside_mode "all" = all cells inside polygon, "distance" = only cells within threshold
#' @return List with profile data frame and celltype contributions
#' @export
calculate_distance_profile <- function(seurat_obj, 
                                       distance_thresholds = seq(10, 200, by = 10),
                                       metric = "all",
                                       group_by = "cellType",
                                       inside_mode = "all") {
  
  if (!inside_mode %in% c("all", "distance")) {
    stop("inside_mode must be either 'all' or 'distance'")
  }
  
  if (!"boundary_distance" %in% colnames(seurat_obj@meta.data)) {
    stop("seurat_obj metadata must contain a 'boundary_distance' column")
  }
  
  cell_names <- colnames(seurat_obj)
  boundary_distance <- seurat_obj$boundary_distance
  n_missing <- sum(is.na(boundary_distance))
  if (n_missing > 0) {
    warning(n_missing, " cell(s) have NA boundary_distance and will be excluded")
  }
  
  results <- data.frame()
  celltype_contributions <- list()
  
  for (dist in distance_thresholds) {
    cat(paste0("Processing distance: ", dist, " μm (inside_mode = '", inside_mode, "')\n"))
    
    if (inside_mode == "all") {
      inside_cells <- cell_names[which(boundary_distance < 0)]
      outside_cells <- cell_names[which(boundary_distance > 0 & 
                                          boundary_distance <= dist)]
    } else if (inside_mode == "distance") {
      inside_cells <- cell_names[which(boundary_distance < 0 & 
                                         abs(boundary_distance) <= dist)]
      outside_cells <- cell_names[which(boundary_distance > 0 & 
                                          boundary_distance <= dist)]
    }
    
    if (length(inside_cells) < 10 || length(outside_cells) < 10) {
      cat(paste0("  Skipping - insufficient cells\n"))
      next
    }
    
    cat(paste0("  Inside cells: ", length(inside_cells), 
               ", Outside cells: ", length(outside_cells), "\n"))
    
    row <- data.frame(distance = dist,
                      n_inside = length(inside_cells),
                      n_outside = length(outside_cells),
                      inside_mode = inside_mode)
    
    if (metric %in% c("composition", "all")) {
      comp_result <- calculate_composition_difference(
        inside_cells, outside_cells, seurat_obj, group_by
      )
      row$composition_difference <- comp_result$total
      
      celltype_contributions[[as.character(dist)]] <- data.frame(
        distance = dist,
        celltype = names(comp_result$by_celltype),
        contribution = comp_result$by_celltype,
        inside_prop = comp_result$inside_prop,
        outside_prop = comp_result$outside_prop,
        diff = comp_result$inside_prop - comp_result$outside_prop
      )
    }
    
    if (metric %in% c("transcriptomic", "all")) {
      row$transcriptomic_distance <- calculate_transcriptomic_distance(
        inside_cells, outside_cells, seurat_obj
      )
    }
    
    if (metric %in% c("de_strength", "all")) {
      row$de_strength <- calculate_de_strength(
        inside_cells, outside_cells, seurat_obj
      )
    }
    
    results <- rbind(results, row)
  }
  
  celltype_df <- do.call(rbind, celltype_contributions)
  
  return(list(
    profile = results,
    celltype_contributions = celltype_df
  ))
}

#' Find Critical Distance
#'
#' Identifies distance at which separation is maximized
#'
#' @param distance_profile Data frame from calculate_distance_profile
#' @param metric_column Column name of metric to analyze
#' @return List with max_distance, max_value, and plateau_distance
#' @export
find_critical_distance <- function(distance_profile, metric_column = "transcriptomic_distance") {
  max_idx <- which.max(distance_profile[[metric_column]])
  max_dist <- distance_profile$distance[max_idx]
  max_val <- distance_profile[[metric_column]][max_idx]
  
  plateau_threshold <- 0.9 * max_val
  plateau_matches <- which(distance_profile[[metric_column]] >= plateau_threshold)
  plateau_dist <- ifelse(length(plateau_matches) > 0, 
                         distance_profile$distance[plateau_matches[1]], 
                         NA)
  
  cat(paste0("\nCritical distances for ", metric_column, ":\n"))
  cat(paste0("  Maximum separation at: ", max_dist, " μm (value: ", round(max_val, 4), ")\n"))
  cat(paste0("  Plateau begins at: ", plateau_dist, " μm\n"))
  
  return(list(max_distance = max_dist, max_value = max_val, plateau_distance = plateau_dist))
}
