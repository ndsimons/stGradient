#' Get Cell Coordinates Keyed by Cell ID
#'
#' Internal helper that returns tissue coordinates matched to Seurat cell IDs.
#' Uses the Seurat v5 `cell` column when present, otherwise the coordinate
#' rownames. Never matches coordinates to cells by position.
#'
#' @param seurat_obj Seurat object.
#' @param fov_name Field of view/image name.
#' @return Data frame with `cell_id`, `x`, `y` and rownames set to cell IDs,
#'   restricted to cells present in both the FOV and `seurat_obj`, in
#'   `colnames(seurat_obj)` order.
#' @keywords internal
get_cell_coords <- function(seurat_obj, fov_name = "fov") {
  coords <- as.data.frame(Seurat::GetTissueCoordinates(seurat_obj, image = fov_name))
  
  if (!all(c("x", "y") %in% colnames(coords))) {
    stop("Coordinates for '", fov_name, "' must contain x/y columns. Found: ",
         paste(colnames(coords), collapse = ", "))
  }
  
  cell_ids <- if ("cell" %in% colnames(coords)) as.character(coords$cell) else rownames(coords)
  out <- data.frame(cell_id = cell_ids, x = coords$x, y = coords$y,
                    stringsAsFactors = FALSE)
  out <- out[!duplicated(out$cell_id), , drop = FALSE]
  rownames(out) <- out$cell_id
  
  shared <- intersect(colnames(seurat_obj), out$cell_id)
  if (length(shared) == 0) {
    stop("No coordinate cell IDs in '", fov_name, "' match colnames(seurat_obj).")
  }
  
  out[shared, , drop = FALSE]
}

# Expand per-cell values to all cells in the object, NA for cells without coordinates.
expand_to_all_cells <- function(values, cell_ids, seurat_obj, what = "values") {
  all_cells <- colnames(seurat_obj)
  n_missing <- length(all_cells) - length(cell_ids)
  if (n_missing > 0) {
    warning(n_missing, " cell(s) have no coordinates in this FOV; ", what,
            " set to NA for those cells.", call. = FALSE)
  }
  out <- values[rep(NA_integer_, length(all_cells))]
  names(out) <- all_cells
  out[cell_ids] <- values
  out
}

# Validate a user-supplied polygon object before spatial operations.
check_polygons <- function(polygons) {
  if (is.null(polygons)) {
    stop("`polygons` is NULL. If they came from add_dbscan_polygons(), retrieve them ",
         "with get_cancer_regions(seurat_obj, fov_name); polygons are no longer stored ",
         "in the FOV boundaries slot.", call. = FALSE)
  }
  if (!inherits(polygons, c("sf", "sfc"))) {
    stop("`polygons` must be an sf or sfc object, not ",
         paste(class(polygons), collapse = "/"), ".", call. = FALSE)
  }
  if (length(sf::st_geometry(polygons)) == 0) {
    stop("`polygons` contains no geometries.", call. = FALSE)
  }
  invisible(polygons)
}

#' Parse SVG Path to Coordinates
#'
#' Internal function to parse plotly SVG paths
#'
#' @param path_string SVG path string from plotly
#' @return Matrix of x,y coordinates
#' @keywords internal
parse_svg_path <- function(path_string) {
  path_clean <- gsub("M|Z", "", path_string)
  
  coords_str <- strsplit(path_clean, "L")[[1]]
  coords_str <- coords_str[coords_str != ""]
  
  coords_list <- lapply(coords_str, function(coord) {
    vals <- as.numeric(strsplit(trimws(coord), ",")[[1]])
    if (length(vals) == 2) return(vals)
    return(NULL)
  })
  
  coords_list <- coords_list[!sapply(coords_list, is.null)]
  
  if (length(coords_list) > 0) {
    coords_matrix <- do.call(rbind, coords_list)
    colnames(coords_matrix) <- c("x", "y")
    
    if (!all(coords_matrix[1,] == coords_matrix[nrow(coords_matrix),])) {
      coords_matrix <- rbind(coords_matrix, coords_matrix[1,])
    }
    
    return(coords_matrix)
  }
  
  return(NULL)
}

#' Add Polygon Membership to Seurat Object
#'
#' Adds metadata column indicating if cells are inside polygons
#'
#' @param seurat_obj Seurat object
#' @param polygons sf polygon object
#' @param fov_name FOV name
#' @return Seurat object with added metadata
#' @export
add_polygon_membership <- function(seurat_obj, polygons, fov_name = "fov") {
  check_polygons(polygons)
  coords <- get_cell_coords(seurat_obj, fov_name)
  coords_sf <- sf::st_as_sf(coords, coords = c("x", "y"))
  
  cells_inside <- sf::st_intersects(coords_sf, polygons, sparse = FALSE)
  is_inside <- apply(cells_inside, 1, any)
  
  seurat_obj$inside_polygon <- expand_to_all_cells(
    is_inside, coords$cell_id, seurat_obj, "inside_polygon"
  )
  
  return(seurat_obj)
}

#' Automatically detect cancer cell regions and create polygons
#'
#' Uses DBSCAN clustering on cancer cell coordinates to identify spatial regions,
#' then creates concave hull polygons around each cluster. Stores the polygons in
#' `seurat_obj@misc$cancer_regions[[fov_name]]` (retrieve them with
#' [get_cancer_regions()]) and adds an `inside_polygon` metadata column.
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param fov_name Name of the field of view (default: "fov")
#' @param cancer_col Metadata column indicating cancer cells (default: "is_cancer")
#' @param eps DBSCAN epsilon parameter - max distance between points (default: 35)
#' @param minPts DBSCAN minimum points per cluster (default: 20)
#' @param concavity Concave hull concavity parameter, lower = tighter fit (default: 2)
#' @param length_threshold Edge length threshold for concave hull (default: 0)
#'
#' @return Seurat object with cancer region polygons stored in `@misc` and
#'   `inside_polygon` metadata
#' @export
#' 
#' @examples
#' \dontrun{
#' seurat_obj <- add_dbscan_polygons(seurat_obj, eps = 50, minPts = 30)
#' }
#'
add_dbscan_polygons <- function(seurat_obj,
                                   fov_name = "fov",
                                   cancer_col = "is_cancer",
                                   eps = 35,
                                   minPts = 20,
                                   concavity = 2,
                                   length_threshold = 0) {
  
  # Check required packages
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }
  if (!requireNamespace("concaveman", quietly = TRUE)) {
    stop("Package 'concaveman' is required. Install with: install.packages('concaveman')")
  }
  
  # Get tissue coordinates
  coords <- get_cell_coords(seurat_obj, fov_name)
  
  # Check if cancer column exists
  if (!cancer_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Column '", cancer_col, "' not found in metadata"))
  }
  
  coords$is_cancer <- seurat_obj@meta.data[coords$cell_id, cancer_col]
  
  # Filter to cancer cells only
  cancer_coords <- coords %>% dplyr::filter(is_cancer == TRUE)
  
  if (nrow(cancer_coords) == 0) {
    stop("No cancer cells found in metadata")
  }
  
  # DBSCAN clustering
  message("Clustering cancer cells with DBSCAN...")
  clusters <- dbscan::dbscan(cancer_coords[, c("x", "y")], 
                             eps = eps, 
                             minPts = minPts)
  cancer_coords$cluster <- clusters$cluster
  
  # Remove noise points (cluster 0)
  cancer_coords <- cancer_coords %>% dplyr::filter(cluster > 0)
  
  if (nrow(cancer_coords) == 0) {
    stop("No clusters found. Try adjusting eps or minPts parameters")
  }
  
  n_clusters <- length(unique(cancer_coords$cluster))
  message(paste0("Found ", n_clusters, " cancer regions"))
  
  # Create polygons for each cluster
  polygon_list <- list()
  for (cluster_id in unique(cancer_coords$cluster)) {
    cluster_coords <- cancer_coords %>% dplyr::filter(cluster == cluster_id)
    cluster_matrix <- as.matrix(cluster_coords[, c("x", "y")])
    
    # Create concave hull
    concave <- concaveman::concaveman(cluster_matrix,
                                      concavity = concavity,
                                      length_threshold = length_threshold)
    
    # Convert to sf polygon
    polygon_list[[cluster_id]] <- sf::st_polygon(list(concave)) %>%
      sf::st_sfc() %>%
      sf::st_sf(region_id = cluster_id)
  }
  
  # Combine all polygons
  auto_polygons <- do.call(rbind, polygon_list)
  
  # Store in @misc: FOV boundaries only accept Segmentation/Centroids objects,
  # and an sf object there breaks subset()
  message("Adding polygons to Seurat object...")
  if (is.null(seurat_obj@misc$cancer_regions)) {
    seurat_obj@misc$cancer_regions <- list()
  }
  seurat_obj@misc$cancer_regions[[fov_name]] <- auto_polygons
  seurat_obj <- add_polygon_membership(seurat_obj, auto_polygons, fov_name = fov_name)
  message(paste0("Successfully added ", n_clusters, " cancer region polygons"))
  
  return(seurat_obj)
}

#' Get Cancer Region Polygons
#'
#' Retrieves the polygons stored by [add_dbscan_polygons()] for a field of view.
#' Also reads polygons from the FOV `boundaries` slot, where objects created by
#' earlier versions of stGradient stored them.
#'
#' @param seurat_obj Seurat object processed with [add_dbscan_polygons()]
#' @param fov_name Name of the field of view (default: "fov")
#' @return sf object of cancer region polygons, or `NULL` with a warning if none
#'   are stored for `fov_name`
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- add_dbscan_polygons(seurat_obj)
#' auto_polygons <- get_cancer_regions(seurat_obj)
#' }
get_cancer_regions <- function(seurat_obj, fov_name = "fov") {
  polygons <- seurat_obj@misc$cancer_regions[[fov_name]]
  if (is.null(polygons) && fov_name %in% names(seurat_obj@images)) {
    legacy <- tryCatch(
      seurat_obj@images[[fov_name]]@boundaries$cancer_regions,
      error = function(e) NULL
    )
    if (inherits(legacy, "sf")) polygons <- legacy
  }
  if (is.null(polygons)) {
    stored <- names(seurat_obj@misc$cancer_regions)
    warning("No cancer region polygons stored for FOV '", fov_name, "'. ",
            if (length(stored) > 0) paste0("Polygons exist for: ", paste(stored, collapse = ", "), ". ")
            else "Run add_dbscan_polygons() first. ",
            call. = FALSE)
  }
  polygons
}

#' Optimize DBSCAN Parameters
#'
#' Performs a grid search over eps and minPts for DBSCAN clustering of cancer
#' cell coordinates. Each combination is scored by average silhouette score and
#' by the fraction of cancer cells labeled noise. Combinations with noise below
#' `noise_pct` are turned into concave hull polygons and scored by coverage
#' (fraction of cancer cells inside a polygon) and purity (fraction of cells
#' inside polygons that are cancer), which are combined into an F-beta score.
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param noise_pct Maximum fraction of cancer cells labeled noise for a
#'   combination to be scored by coverage and purity (default: 0.15)
#' @param b Beta for the F-beta score. Values below 1 weight purity more than
#'   coverage; values above 1 weight coverage more (default: 0.5)
#' @param plot Print a heatmap of coverage and purity for the scored
#'   combinations (default: TRUE)
#' @param fov_name Name of the field of view (default: "fov")
#' @param cancer_col Metadata column indicating cancer cells (default: "is_cancer")
#' @param eps_range Numeric vector of eps values to test (default: seq(10, 100, by = 5))
#' @param minPts_range Numeric vector of minPts values to test (default: c(10, 15, 20, 30, 50))
#' @param concavity Concave hull concavity parameter used to build polygons (default: 2)
#' @param length_threshold Edge length threshold for concave hull (default: 0)
#' @param max_sample Maximum number of non-noise points to sample for silhouette
#'   calculation to control memory usage (default: 46000)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A data frame with columns eps, minPts, sil_score, n_clusters,
#'   noise_pct, n_regions, cancer_covered, purity and fbeta, sorted by
#'   descending fbeta (then silhouette score). The best combination is in the
#'   first row. Combinations above `noise_pct` have `NA` coverage, purity and
#'   fbeta. When `plot = TRUE`, the heatmap is also stored in `attr(result, "plot")`.
#' @export
#'
#' @examples
#' \dontrun{
#' results <- optimize_dbscan_params(seurat_obj, noise_pct = 0.15, b = 0.5, plot = TRUE)
#' head(results)
#'
#' # Use best parameters
#' best <- results[1, ]
#' seurat_obj <- add_dbscan_polygons(seurat_obj, eps = best$eps, minPts = best$minPts)
#' }
#'
optimize_dbscan_params <- function(seurat_obj,
                                   noise_pct = 0.15,
                                   b = 0.5,
                                   plot = TRUE,
                                   fov_name = "fov",
                                   cancer_col = "is_cancer",
                                   eps_range = seq(10, 100, by = 5),
                                   minPts_range = c(10, 15, 20, 30, 50),
                                   concavity = 2,
                                   length_threshold = 0,
                                   max_sample = 46000,
                                   verbose = TRUE) {

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required. Install with: install.packages('cluster')")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }
  if (!is.numeric(noise_pct) || length(noise_pct) != 1 || noise_pct < 0 || noise_pct > 1) {
    stop("`noise_pct` must be a single number between 0 and 1")
  }
  if (!is.numeric(b) || length(b) != 1 || b <= 0) {
    stop("`b` must be a single positive number")
  }

  coords <- get_cell_coords(seurat_obj, fov_name)

  if (!cancer_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Column '", cancer_col, "' not found in metadata"))
  }

  coords$is_cancer <- as.logical(seurat_obj@meta.data[coords$cell_id, cancer_col])
  coords$is_cancer[is.na(coords$is_cancer)] <- FALSE
  cancer_coords <- coords[coords$is_cancer, , drop = FALSE]

  if (nrow(cancer_coords) == 0) {
    stop("No cancer cells found in metadata")
  }

  xy <- cancer_coords[, c("x", "y")]
  all_pts <- sf::st_as_sf(coords, coords = c("x", "y"))

  results <- expand.grid(eps = eps_range, minPts = minPts_range)
  results$sil_score <- NA_real_
  results$n_clusters <- NA_integer_
  results$noise_pct <- NA_real_
  results$n_regions <- NA_integer_
  results$cancer_covered <- NA_real_
  results$purity <- NA_real_

  n_combos <- nrow(results)
  if (verbose) message(sprintf("Testing %d parameter combinations...", n_combos))

  for (i in seq_len(n_combos)) {
    if (verbose) {
      message(sprintf("  %d/%d  eps=%.0f  minPts=%d",
                      i, n_combos, results$eps[i], results$minPts[i]))
    }

    cl <- dbscan::dbscan(xy, eps = results$eps[i], minPts = results$minPts[i])
    results$n_clusters[i] <- max(cl$cluster)
    results$noise_pct[i] <- mean(cl$cluster == 0)

    # Need at least 2 clusters and more points than clusters for silhouette
    if (max(cl$cluster) >= 2 && sum(cl$cluster > 0) > max(cl$cluster)) {
      non_noise <- which(cl$cluster > 0)

      # Subsample to limit memory for distance matrix
      if (length(non_noise) > max_sample) {
        non_noise <- sample(non_noise, max_sample)
      }

      sil <- cluster::silhouette(cl$cluster[non_noise], stats::dist(xy[non_noise, ]))
      results$sil_score[i] <- mean(sil[, 3])
    }

    # Coverage and purity of the resulting polygons for low-noise combinations
    if (max(cl$cluster) >= 1 && results$noise_pct[i] < noise_pct) {
      polys <- lapply(sort(unique(cl$cluster[cl$cluster > 0])), function(k) {
        hull <- concaveman::concaveman(as.matrix(xy[cl$cluster == k, ]),
                                       concavity = concavity,
                                       length_threshold = length_threshold)
        sf::st_polygon(list(hull))
      })
      region <- sf::st_union(sf::st_make_valid(sf::st_sfc(polys)))
      inside <- lengths(sf::st_intersects(all_pts, region)) > 0
      results$n_regions[i] <- length(polys)
      results$cancer_covered[i] <- mean(inside[coords$is_cancer])
      results$purity[i] <- if (any(inside)) mean(coords$is_cancer[inside]) else NA_real_
    }
  }

  results$fbeta <- with(results, (1 + b^2) * purity * cancer_covered /
                          (b^2 * purity + cancer_covered))

  results <- results[order(-results$fbeta, -results$sil_score, na.last = TRUE), ]
  rownames(results) <- NULL

  scored <- results[!is.na(results$cancer_covered), , drop = FALSE]

  if (verbose) {
    if (nrow(scored) > 0) {
      best <- results[1, ]
      message(sprintf(
        "Best: eps=%.0f, minPts=%d (fbeta=%.3f, coverage=%.1f%%, purity=%.1f%%, %.1f%% noise)",
        best$eps, best$minPts, best$fbeta, best$cancer_covered * 100,
        best$purity * 100, best$noise_pct * 100))
    } else {
      message(sprintf(
        "No combinations had noise below %.0f%%. Try a larger noise_pct or different parameter ranges.",
        noise_pct * 100))
    }
  }

  if (isTRUE(plot) && nrow(scored) > 0) {
    long <- tidyr::pivot_longer(scored, c("cancer_covered", "purity"),
                                names_to = "metric", values_to = "value")
    p <- ggplot2::ggplot(long, ggplot2::aes(factor(eps), factor(minPts), fill = value)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)), size = 2.8) +
      ggplot2::facet_wrap(~ metric, ncol = 1) +
      ggplot2::scale_fill_viridis_c(labels = function(x) sprintf("%.0f%%", 100 * x)) +
      ggplot2::labs(x = "eps", y = "minPts", fill = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    print(p)
    attr(results, "plot") <- p
  }

  results
}

#' Plot DBSCAN-based cancer region polygons
#'
#' Uses DBSCAN clustering on cancer cell coordinates to identify spatial regions,
#' creates concave hull polygons around each cluster, and optionally visualizes them.
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param fov Name of the field of view (default: "fov")
#' @param eps DBSCAN epsilon parameter - max distance between points (default: 35)
#' @param minPts DBSCAN minimum points per cluster (default: 20)
#' @param concavity Concave hull concavity parameter, lower = tighter fit (default: 2)
#' @param length_threshold Edge length threshold for concave hull (default: 0)
#' @param sample_name Optional sample name prefixed to the plot titles
#' @param plot Whether to generate visualization plots (default: TRUE)
#' @param ncol Number of plot columns: 2 places the plots side by side, 1 stacks
#'   them vertically (default: 2)
#'
#' @return sf object containing polygon geometries for each cancer region
#' @export
#' 
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom sf st_as_sf st_polygon st_sfc st_sf
#' @importFrom dbscan dbscan
#' @importFrom concaveman concaveman
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_sf aes scale_color_manual scale_x_reverse scale_y_reverse theme_minimal theme element_rect element_blank element_text ggtitle
#' @importFrom cowplot plot_grid
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' polygons <- plotDBSCANpolygons(seurat_obj, sample_name = "Patient1")
#' 
#' # Custom parameters without plotting
#' polygons <- plotDBSCANpolygons(
#'   seurat_obj, 
#'   eps = 50, 
#'   minPts = 30,
#'   sample_name = "Patient1",
#'   plot = FALSE
#' )
#' }
#'
plot_dbscan_polygons <- function(seurat_obj, 
                                fov = "fov", 
                                eps = 35, 
                                minPts = 20, 
                                concavity = 2,
                                length_threshold = 0,
                                sample_name = NULL,
                                plot = TRUE,
                                ncol = 2) {
  if (!ncol %in% c(1, 2)) {
    stop("`ncol` must be 1 or 2")
  }
  
  # Check required packages
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }
  if (!requireNamespace("concaveman", quietly = TRUE)) {
    stop("Package 'concaveman' is required. Install with: install.packages('concaveman')")
  }
  
  # Get coordinates and cancer annotation
  coords <- get_cell_coords(seurat_obj, fov)
  coords$is_cancer <- seurat_obj@meta.data[coords$cell_id, "is_cancer"]
  
  # Convert to sf object
  coords_sf <- sf::st_as_sf(coords, coords = c("x", "y"))
  
  # Filter cancer cells
  cancer_coords <- coords %>% dplyr::filter(is_cancer == TRUE)
  
  if (nrow(cancer_coords) == 0) {
    stop("No cancer cells found in seurat_obj$is_cancer")
  }
  
  # DBSCAN clustering
  clusters <- dbscan::dbscan(cancer_coords[, c("x", "y")], eps = eps, minPts = minPts)
  cancer_coords$cluster <- clusters$cluster
  
  # Remove noise points (cluster 0)
  cancer_coords <- cancer_coords %>% dplyr::filter(cluster > 0)
  
  if (nrow(cancer_coords) == 0) {
    warning("No clusters found. Try adjusting eps or minPts parameters")
    return(NULL)
  }
  
  # Create polygons for each cluster
  polygon_list <- list()
  for (cluster_id in unique(cancer_coords$cluster)) {
    cluster_coords <- cancer_coords %>% dplyr::filter(cluster == cluster_id)
    cluster_matrix <- as.matrix(cluster_coords[, c("x", "y")])
    
    # Create concave hull
    concave <- concaveman::concaveman(cluster_matrix, 
                                      concavity = concavity,
                                      length_threshold = length_threshold)
    
    # Convert to sf polygon
    polygon_list[[cluster_id]] <- sf::st_polygon(list(concave)) %>%
      sf::st_sfc() %>%
      sf::st_sf(region_id = cluster_id)
  }
  
  # Combine all polygons
  auto_polygons <- do.call(rbind, polygon_list)
  
  # Visualize if requested
  if (plot) {
    title_prefix <- if (is.null(sample_name)) "" else paste0(sample_name, ": ")
    
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = coords_sf, ggplot2::aes(color = is_cancer), size = 0.005, alpha = 0.5) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "lightgray")) +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"),
                     plot.background = ggplot2::element_rect(fill = "black"),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(color = "white")) +
      ggplot2::ggtitle(paste0(title_prefix, "cells of interest"))
    
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = coords_sf, ggplot2::aes(color = is_cancer), size = 0.005, alpha = 0.5) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "lightgray")) +
      ggplot2::geom_sf(data = auto_polygons, 
                       fill = NA,
                       color = "blue",
                       linewidth = 1) +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"),
                     plot.background = ggplot2::element_rect(fill = "black"),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(color = "white")) +
      ggplot2::ggtitle(paste0(title_prefix, "regions of interest (ROIs)"))
    
    print(cowplot::plot_grid(p1, p2, ncol = ncol))
  }
  
  # Return the polygons
  return(auto_polygons)
}

#' Load Xenium Seurat Object with Stable Cell IDs
#'
#' Reads a Seurat object from RDS, enforces unique cell IDs, optionally computes
#' cancer and ROI metadata columns, and stores coordinates keyed by cell ID.
#'
#' @param rds_path Path to serialized Seurat object (.rds)
#' @param image Optional image/FOV name to use for coordinates
#' @param celltype_col Metadata column containing cell type annotations
#' @param cancer_labels Optional character vector of labels to mark as cancer
#' @param roi_cells Optional character vector of cell IDs to mark as inside ROI
#' @param roi_col Metadata column name for ROI membership
#' @param store_coords Logical; store processed coordinate data in misc
#' @param coords_name Name used for stored coordinate data in misc
#' @return Seurat object with standardized metadata and optional stored coordinates
#' @export
load_xenium_seurat <- function(
    rds_path,
    image = NULL,
    celltype_col,
    cancer_labels = NULL,
    roi_cells = NULL,
    roi_col = "inside_polygon",
    store_coords = TRUE,
    coords_name = "tissue_coords"
) {
  stopifnot(file.exists(rds_path))
  obj <- readRDS(rds_path)

  if (is.null(colnames(obj)) || length(colnames(obj)) == 0) {
    stop("Seurat object has no colnames (cell IDs). Cannot proceed.")
  }
  if (anyDuplicated(colnames(obj)) > 0) {
    original <- colnames(obj)
    colnames(obj) <- make.unique(original)
    obj$cell_id_original <- original
  }

  obj$cell_id <- colnames(obj)

  if (!is.null(cancer_labels)) {
    if (celltype_col %in% colnames(obj@meta.data)) {
      obj$is_cancer <- obj[[celltype_col, drop = TRUE]] %in% cancer_labels
    } else {
      warning(
        "cancer_labels provided, but celltype_col='", celltype_col,
        "' not found in meta.data. Skipping is_cancer."
      )
    }
  }

  if (!is.null(roi_cells)) {
    roi_cells <- intersect(roi_cells, colnames(obj))
    obj[[roi_col]] <- colnames(obj) %in% roi_cells
  }

  imgs <- names(obj@images)
  if (length(imgs) > 0) {
    if (is.null(image)) {
      image <- imgs[1]
    } else if (!image %in% imgs) {
      stop("Requested image not found. Available images: ", paste(imgs, collapse = ", "))
    }

    df <- get_cell_coords(obj, image)

    if (roi_col %in% colnames(obj@meta.data)) {
      df[[roi_col]] <- obj@meta.data[df$cell_id, roi_col]
    }
    if ("is_cancer" %in% colnames(obj@meta.data)) {
      df$is_cancer <- obj@meta.data[df$cell_id, "is_cancer"]
    }
    if (celltype_col %in% colnames(obj@meta.data)) {
      df[[celltype_col]] <- obj@meta.data[df$cell_id, celltype_col]
    }

    obj@misc$xenium_image_name <- image
    if (store_coords) {
      obj@misc[[coords_name]] <- df
    }
  } else {
    warning("No images/FOV found in obj@images; skipping coordinate extraction.")
  }

  return(obj)
}
