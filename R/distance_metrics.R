#' Calculate Boundary Distances
#'
#' Calculates signed distance from each cell to polygon boundary
#'
#' @param seurat_obj Seurat object with spatial data
#' @param polygons sf polygon object
#' @param fov_name Name of FOV in Seurat object
#' @return Named numeric vector of signed distances (negative = inside), in
#'   `colnames(seurat_obj)` order; `NA` for cells without coordinates in the FOV.
#' @export
calculate_boundary_distances <- function(seurat_obj, polygons, fov_name = "fov") {
  coords <- get_cell_coords(seurat_obj, fov_name)
  coords_sf <- sf::st_as_sf(coords, coords = c("x", "y"))
  
  polygon_boundaries <- sf::st_cast(polygons, "MULTILINESTRING")
  
  distances <- sf::st_distance(coords_sf, polygon_boundaries)
  min_distances <- apply(distances, 1, min)
  
  cells_inside <- sf::st_intersects(coords_sf, polygons, sparse = FALSE)
  is_inside <- apply(cells_inside, 1, any)
  
  signed_distances <- ifelse(is_inside, -min_distances, min_distances)
  
  return(expand_to_all_cells(signed_distances, coords$cell_id, seurat_obj, "boundary distances"))
}

#' Calculate Transcriptomic Distance
#'
#' Compares average gene expression profiles between two cell groups
#'
#' @param inside_cells Character vector of cell IDs inside region
#' @param outside_cells Character vector of cell IDs outside region
#' @param seurat_obj Seurat object
#' @return Numeric value (1 - correlation)
#' @export
calculate_transcriptomic_distance <- function(inside_cells, outside_cells, seurat_obj) {
  # Prefer Seurat v5 layer API and fall back for older versions.
  expr_data <- tryCatch(
    Seurat::GetAssayData(seurat_obj, layer = "data"),
    error = function(e) Seurat::GetAssayData(seurat_obj, slot = "data")
  )
  
  # Subset to cells of interest - force to matrix
  inside_expr <- as.matrix(expr_data[, inside_cells, drop = FALSE])
  outside_expr <- as.matrix(expr_data[, outside_cells, drop = FALSE])
  
  # Calculate mean profiles
  inside_profile <- rowMeans(inside_expr)
  outside_profile <- rowMeans(outside_expr)
  
  # Calculate distance as 1 - correlation
  dist <- 1 - stats::cor(inside_profile, outside_profile)
  return(dist)
}

#' Calculate Differential Expression Strength
#'
#' Calculates mean absolute log2 fold change of top DE genes
#'
#' @param inside_cells Character vector of cell IDs inside region
#' @param outside_cells Character vector of cell IDs outside region
#' @param seurat_obj Seurat object
#' @return Numeric value (mean |log2FC|)
#' @export
calculate_de_strength <- function(inside_cells, outside_cells, seurat_obj) {
  if (length(intersect(inside_cells, outside_cells)) > 0) {
    stop("`inside_cells` and `outside_cells` must not overlap.")
  }
  
  # Cells in neither group are labelled "excluded" so they don't enter the comparison
  cell_ids <- colnames(seurat_obj)
  temp_ident <- rep("excluded", length(cell_ids))
  temp_ident[cell_ids %in% inside_cells] <- "inside"
  temp_ident[cell_ids %in% outside_cells] <- "outside"
  seurat_obj$temp_group <- temp_ident
  
  Seurat::Idents(seurat_obj) <- "temp_group"
  de_genes <- Seurat::FindMarkers(seurat_obj, 
                                   ident.1 = "inside", 
                                   ident.2 = "outside",
                                   logfc.threshold = 0,
                                   min.pct = 0.1,
                                   verbose = FALSE)
  
  if (nrow(de_genes) > 0) {
    top_genes <- utils::head(de_genes[order(abs(de_genes$avg_log2FC), decreasing = TRUE), ], 100)
    strength <- mean(abs(top_genes$avg_log2FC))
  } else {
    strength <- 0
  }
  
  return(strength)
}

#' Calculate Cell Type Composition Difference
#'
#' Compares cell type distributions between two groups
#'
#' @param inside_cells Character vector of cell IDs inside region
#' @param outside_cells Character vector of cell IDs outside region
#' @param seurat_obj Seurat object
#' @param group_by Metadata column with cell type annotations
#' @return List with total difference and per-celltype breakdown
#' @export
calculate_composition_difference <- function(inside_cells, outside_cells, seurat_obj, group_by = "cellType") {
  inside_types <- table(seurat_obj[[group_by]][inside_cells, ])
  outside_types <- table(seurat_obj[[group_by]][outside_cells, ])
  
  all_types <- unique(c(names(inside_types), names(outside_types)))
  inside_vec <- sapply(all_types, function(x) ifelse(x %in% names(inside_types), inside_types[x], 0))
  outside_vec <- sapply(all_types, function(x) ifelse(x %in% names(outside_types), outside_types[x], 0))
  
  inside_prop <- inside_vec / sum(inside_vec)
  outside_prop <- outside_vec / sum(outside_vec)
  
  sq_diff <- (inside_prop - outside_prop)^2
  names(sq_diff) <- all_types
  
  total_diff <- sum(sq_diff)
  
  return(list(
    total = total_diff,
    by_celltype = sq_diff,
    inside_prop = inside_prop,
    outside_prop = outside_prop
  ))
}
