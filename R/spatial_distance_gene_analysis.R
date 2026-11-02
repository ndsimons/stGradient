#' Compute Signed Distance From Tumor Boundary
#'
#' Computes per-cell signed distance to an intratumoral boundary polygon built
#' from intratumoral cells using a concave hull.
#'
#' @param seurat_obj Seurat object.
#' @param polygon_col Metadata column indicating inside/outside polygon.
#' @param intratumoral_label Value in `polygon_col` that marks intratumoral cells.
#' @param fov_name Field of view/image name for coordinate extraction.
#' @param concavity Concavity parameter passed to `concaveman`.
#' @return Seurat object with `dist_to_boundary` and `abs_dist_to_boundary` metadata.
#' @export
compute_distance_from_boundary <- function(seurat_obj,
                                           polygon_col = "inside_polygon",
                                           intratumoral_label = TRUE,
                                           fov_name = "fov",
                                           concavity = 2) {
  coords <- get_cell_coords(seurat_obj, fov_name)

  meta <- seurat_obj@meta.data[coords$cell_id, , drop = FALSE]
  is_intratumoral <- as.character(meta[[polygon_col]]) == as.character(intratumoral_label)

  if (sum(is_intratumoral, na.rm = TRUE) < 10) {
    stop("Too few intratumoral cells to build a boundary polygon.")
  }

  intra_coords <- coords[which(is_intratumoral), c("x", "y"), drop = FALSE]
  intra_sf <- sf::st_as_sf(intra_coords, coords = c("x", "y"))
  boundary_polygon <- concaveman::concaveman(intra_sf, concavity = concavity)

  all_points <- sf::st_as_sf(coords, coords = c("x", "y"))
  boundary_line <- sf::st_boundary(boundary_polygon)
  distances <- as.numeric(sf::st_distance(all_points, boundary_line))
  inside <- sf::st_within(all_points, boundary_polygon, sparse = FALSE)[, 1]
  signed_distance <- ifelse(inside, -distances, distances)
  signed_distance <- expand_to_all_cells(signed_distance, coords$cell_id, seurat_obj,
                                         "dist_to_boundary")

  seurat_obj$dist_to_boundary <- signed_distance
  seurat_obj$abs_dist_to_boundary <- abs(signed_distance)
  seurat_obj
}

# Fetch the normalized "data" layer, running NormalizeData() if it is empty.
get_data_layer <- function(obj, assay) {
  if (!assay %in% names(obj@assays)) {
    stop("Assay '", assay, "' not found. Available assays: ",
         paste(names(obj@assays), collapse = ", "), call. = FALSE)
  }
  fetch <- function(o) suppressWarnings(tryCatch(
    Seurat::GetAssayData(o, assay = assay, layer = "data"),
    error = function(e) Seurat::GetAssayData(o, assay = assay, slot = "data")
  ))
  expr_mat <- fetch(obj)
  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    message("The 'data' layer of assay '", assay, "' is empty; running ",
            "Seurat::NormalizeData() with default settings.")
    obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)
    expr_mat <- fetch(obj)
  }
  expr_mat
}

# Warn why gene filtering left too few genes to test.
warn_too_few_genes <- function(n_expressing, min_cells_expressing, n_cells, assay) {
  max_expr <- if (length(n_expressing) > 0) max(n_expressing) else 0
  warning(sprintf(paste0(
    "Only %d gene(s) expressed in >= %d of %d cells in the '%s' assay 'data' layer ",
    "(%d genes total; max cells expressing any gene: %d); returning NULL. ",
    "Try lowering min_cells_expressing."),
    sum(n_expressing >= min_cells_expressing), min_cells_expressing, n_cells, assay,
    length(n_expressing), max_expr
  ), call. = FALSE)
}

#' Find Distance-Associated Genes
#'
#' Fits per-gene linear models and Spearman correlations against distance from the
#' tumor boundary.
#'
#' @param seurat_obj Seurat object with `dist_to_boundary` metadata.
#' @param assay Assay name. If its `data` layer is empty, `Seurat::NormalizeData()`
#'   is run with default settings on the analyzed cells.
#' @param min_cells_expressing Minimum expressing cells for a gene to be tested.
#' @param n_distance_bins Number of bins for trend summaries.
#' @param focus Region to analyze: `peritumoral`, `intratumoral`, or `both`.
#' @param max_dist_quantile Distance quantile used to trim extreme outliers.
#' @param p_adj_threshold BH adjusted p-value threshold for significant genes.
#' @return List containing model results, binned expression, and metadata.
#' @export
find_distance_associated_genes <- function(seurat_obj,
                                           assay = "Xenium",
                                           min_cells_expressing = 10,
                                           n_distance_bins = 10,
                                           focus = c("peritumoral", "intratumoral", "both"),
                                           max_dist_quantile = 0.95,
                                           p_adj_threshold = 0.05) {
  focus <- match.arg(focus)

  if (!"dist_to_boundary" %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column 'dist_to_boundary' not found. Run compute_distance_from_boundary() first.")
  }

  dist_all <- seurat_obj$dist_to_boundary
  has_dist <- !is.na(dist_all)
  if (focus == "peritumoral") {
    cells_use <- colnames(seurat_obj)[has_dist & dist_all > 0]
  } else if (focus == "intratumoral") {
    cells_use <- colnames(seurat_obj)[has_dist & dist_all <= 0]
  } else {
    cells_use <- colnames(seurat_obj)[has_dist]
  }

  if (length(cells_use) < 20) {
    warning(sprintf(
      "Only %d %s cells (need >= 20); returning NULL. dist_to_boundary: %d > 0, %d <= 0, %d NA.",
      length(cells_use), focus, sum(dist_all > 0, na.rm = TRUE),
      sum(dist_all <= 0, na.rm = TRUE), sum(!has_dist)
    ), call. = FALSE)
    return(NULL)
  }

  sub_obj <- subset(seurat_obj, cells = cells_use)

  if (focus == "intratumoral") {
    dist_values <- abs(sub_obj$dist_to_boundary)
  } else {
    dist_values <- sub_obj$dist_to_boundary
  }

  max_dist <- stats::quantile(dist_values, max_dist_quantile, na.rm = TRUE)
  keep <- dist_values <= max_dist
  sub_obj <- subset(sub_obj, cells = colnames(sub_obj)[keep])
  dist_values <- dist_values[keep]

  expr_mat <- get_data_layer(sub_obj, assay)

  n_expressing <- Matrix::rowSums(expr_mat > 0)
  genes_use <- names(n_expressing[n_expressing >= min_cells_expressing])

  if (length(genes_use) < 2) {
    warn_too_few_genes(n_expressing, min_cells_expressing, ncol(expr_mat), assay)
    return(NULL)
  }

  expr_mat <- expr_mat[genes_use, , drop = FALSE]

  results_lm <- data.frame(
    gene = genes_use,
    coef = NA_real_,
    se = NA_real_,
    t_stat = NA_real_,
    p_value = NA_real_,
    r_squared = NA_real_,
    mean_expr = NA_real_,
    pct_expressing = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes_use)) {
    gene <- genes_use[i]
    y <- as.numeric(expr_mat[gene, ])
    fit <- tryCatch(summary(stats::lm(y ~ dist_values)), error = function(e) NULL)
    if (!is.null(fit)) {
      results_lm$coef[i] <- fit$coefficients[2, 1]
      results_lm$se[i] <- fit$coefficients[2, 2]
      results_lm$t_stat[i] <- fit$coefficients[2, 3]
      results_lm$p_value[i] <- fit$coefficients[2, 4]
      results_lm$r_squared[i] <- fit$r.squared
    }
    results_lm$mean_expr[i] <- mean(y)
    results_lm$pct_expressing[i] <- mean(y > 0) * 100
  }

  results_lm$p_adj <- stats::p.adjust(results_lm$p_value, method = "BH")
  results_lm$direction <- ifelse(results_lm$coef > 0, "increases_with_distance", "decreases_with_distance")
  results_lm <- results_lm %>% dplyr::arrange(p_adj)

  bin_breaks <- unique(stats::quantile(
    dist_values,
    probs = seq(0, 1, length.out = n_distance_bins + 1),
    na.rm = TRUE
  ))
  dist_bins <- cut(dist_values, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)
  bin_midpoints <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2

  sig_genes <- results_lm$gene[results_lm$p_adj < p_adj_threshold]
  if (length(sig_genes) > 0) {
    top_genes <- utils::head(sig_genes, min(50, length(sig_genes)))
    binned_expr <- do.call(rbind, lapply(top_genes, function(gene) {
      y <- as.numeric(expr_mat[gene, ])
      bvals <- tapply(y, dist_bins, mean)
      data.frame(
        gene = gene,
        bin = seq_along(bvals),
        bin_midpoint = bin_midpoints[seq_along(bvals)],
        mean_expr = as.numeric(bvals),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    binned_expr <- NULL
  }

  spearman_results <- data.frame(
    gene = genes_use,
    rho = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes_use)) {
    gene <- genes_use[i]
    y <- as.numeric(expr_mat[gene, ])
    cor_test <- tryCatch(
      stats::cor.test(dist_values, y, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(cor_test)) {
      spearman_results$rho[i] <- unname(cor_test$estimate)
      spearman_results$p_value[i] <- cor_test$p.value
    }
  }
  spearman_results$p_adj <- stats::p.adjust(spearman_results$p_value, method = "BH")

  results <- merge(results_lm, spearman_results, by = "gene", suffixes = c("_lm", "_spearman"))
  results <- results %>% dplyr::arrange(p_adj_lm)

  list(
    results = results,
    binned_expr = binned_expr,
    dist_values = dist_values,
    focus = focus,
    n_distance_bins = n_distance_bins
  )
}

#' Plot Distance-Associated Gene Results
#'
#' @param result_list Output list from `find_distance_associated_genes`.
#' @param top_n Maximum number of significant genes shown in trend facets.
#' @param sample_label Optional label for plot titles.
#' @return List with `volcano` and `trends` ggplot objects.
#' @export
plot_distance_genes <- function(result_list, top_n = 12, sample_label = "") {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")
  }

  results <- result_list$results
  binned <- result_list$binned_expr
  focus <- result_list$focus

  title_suffix <- ifelse(
    nchar(sample_label) > 0,
    sprintf(" - %s (%s)", sample_label, focus),
    sprintf(" (%s)", focus)
  )

  p_volcano <- ggplot2::ggplot(results, ggplot2::aes(x = coef, y = -log10(p_adj_lm))) +
    ggplot2::geom_point(ggplot2::aes(color = direction), alpha = 0.5, size = 1) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(
      data = results %>%
        dplyr::filter(p_adj_lm < 0.05) %>%
        dplyr::group_by(direction) %>%
        dplyr::slice_min(p_adj_lm, n = 5),
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    ggplot2::scale_color_manual(values = c(
      increases_with_distance = "#E74C3C",
      decreases_with_distance = "#3498DB"
    )) +
    ggplot2::labs(
      title = paste0("Distance-associated genes", title_suffix),
      x = "Coefficient (expression ~ distance)",
      y = "-log10(adjusted p-value)"
    ) +
    ggplot2::theme_minimal()

  p_trends <- NULL
  if (!is.null(binned)) {
    n_up <- ceiling(top_n / 2)
    n_down <- top_n - n_up
    top_up <- results %>%
      dplyr::filter(p_adj_lm < 0.05, coef > 0) %>%
      dplyr::slice_min(p_adj_lm, n = n_up) %>%
      dplyr::pull(gene)
    top_down <- results %>%
      dplyr::filter(p_adj_lm < 0.05, coef < 0) %>%
      dplyr::slice_min(p_adj_lm, n = n_down) %>%
      dplyr::pull(gene)
    plot_data <- binned %>% dplyr::filter(gene %in% c(top_up, top_down))

    if (nrow(plot_data) > 0) {
      p_trends <- ggplot2::ggplot(plot_data, ggplot2::aes(x = bin_midpoint, y = mean_expr, color = gene)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 4) +
        ggplot2::labs(
          title = paste0("Expression trends", title_suffix),
          x = "Distance from boundary",
          y = "Mean expression"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    }
  }

  list(volcano = p_volcano, trends = p_trends)
}

#' Define Intratumoral T Cell Hotspots
#'
#' Uses DBSCAN on intratumoral CD4+/CD8+ T cells and returns hotspot polygons,
#' per-cell hotspot assignment, and per-cell distance to nearest hotspot.
#'
#' @param seurat_obj Seurat object.
#' @param celltype_col Metadata cell type column.
#' @param tcell_labels Cell type labels considered T cells.
#' @param polygon_col Metadata polygon membership column.
#' @param intratumoral_label Intratumoral value in `polygon_col`.
#' @param eps DBSCAN neighborhood radius.
#' @param min_pts DBSCAN minimum points.
#' @param min_hotspot_size Minimum cells to keep a hotspot cluster.
#' @param fov_name Field of view/image name for coordinate extraction.
#' @return List with hotspot polygons, per-cell metadata, and summary fields.
#' @export
define_tcell_hotspots <- function(seurat_obj,
                                  celltype_col = "singleR.predicted.id.brcaAtlas",
                                  tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
                                  polygon_col = "inside_polygon",
                                  intratumoral_label = TRUE,
                                  eps = 50,
                                  min_pts = 10,
                                  min_hotspot_size = 5,
                                  fov_name = "fov") {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }

  coords <- get_cell_coords(seurat_obj, fov_name)
  meta <- seurat_obj@meta.data[coords$cell_id, , drop = FALSE]
  coords <- coords[, c("x", "y"), drop = FALSE]

  is_intra <- as.character(meta[[polygon_col]]) == as.character(intratumoral_label)
  is_tcell <- as.character(meta[[celltype_col]]) %in% tcell_labels
  keep_mask <- is_intra & is_tcell & !is.na(is_intra) & !is.na(is_tcell)

  if (sum(keep_mask) < min_pts * 2) {
    return(NULL)
  }

  tcell_coords <- coords[keep_mask, , drop = FALSE]
  db <- dbscan::dbscan(as.matrix(tcell_coords), eps = eps, minPts = min_pts)
  tcell_coords$hotspot_cluster <- db$cluster

  cluster_sizes <- table(tcell_coords$hotspot_cluster)
  valid_clusters <- as.integer(names(cluster_sizes[
    names(cluster_sizes) != "0" & cluster_sizes >= min_hotspot_size
  ]))

  if (length(valid_clusters) == 0) {
    return(NULL)
  }

  hotspot_polys <- lapply(valid_clusters, function(cl) {
    pts <- tcell_coords[tcell_coords$hotspot_cluster == cl, c("x", "y"), drop = FALSE]
    if (nrow(pts) < 3) {
      return(NULL)
    }
    hull_idx <- grDevices::chull(pts)
    hull_pts <- pts[c(hull_idx, hull_idx[1]), , drop = FALSE]
    sf::st_polygon(list(as.matrix(hull_pts)))
  })
  hotspot_polys <- Filter(Negate(is.null), hotspot_polys)

  hotspot_sf <- sf::st_sf(
    hotspot_id = seq_along(hotspot_polys),
    geometry = sf::st_sfc(hotspot_polys)
  )

  all_pts_sf <- sf::st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]), coords = c("x", "y"))
  hotspot_union <- sf::st_union(hotspot_sf)
  dist_to_hs <- as.numeric(sf::st_distance(all_pts_sf, hotspot_union))
  inside_hs <- sf::st_within(all_pts_sf, hotspot_union, sparse = FALSE)[, 1]
  signed_dist <- ifelse(inside_hs, -dist_to_hs, dist_to_hs)

  cluster_lookup <- rep(NA_integer_, ncol(seurat_obj))
  names(cluster_lookup) <- colnames(seurat_obj)
  cluster_lookup[rownames(tcell_coords)] <- tcell_coords$hotspot_cluster

  hotspot_assigned <- stats::setNames(rep(FALSE, ncol(seurat_obj)), colnames(seurat_obj))
  hotspot_assigned[rownames(coords)] <- keep_mask

  tcell_meta <- data.frame(
    cell_id = colnames(seurat_obj),
    hotspot_assigned = unname(hotspot_assigned),
    hotspot_cluster = unname(cluster_lookup),
    dist_to_hotspot = unname(expand_to_all_cells(dist_to_hs, rownames(coords), seurat_obj,
                                                 "hotspot distances")),
    signed_dist_hotspot = unname(suppressWarnings(
      expand_to_all_cells(signed_dist, rownames(coords), seurat_obj)
    )),
    stringsAsFactors = FALSE
  )

  list(
    hotspot_sf = hotspot_sf,
    tcell_meta = tcell_meta,
    n_hotspots = length(valid_clusters),
    tcell_coords = tcell_coords
  )
}

#' Find Genes Associated With Distance From T Cell Hotspots
#'
#' @param seurat_obj Seurat object.
#' @param hotspot_result Output from `define_tcell_hotspots`.
#' @param assay Assay name. If its `data` layer is empty, `Seurat::NormalizeData()`
#'   is run with default settings on the analyzed cells.
#' @param focus Region to analyze relative to hotspot: `perihotspot`,
#'   `intrahotspot`, or `both`.
#' @param min_cells_expressing Minimum expressing cells for a gene to be tested.
#' @param n_distance_bins Number of bins for trend summaries.
#' @param max_dist_quantile Distance quantile used to trim extreme outliers.
#' @param p_adj_threshold BH adjusted p-value threshold.
#' @param exclude_tcell_labels Optional labels to exclude from analysis.
#' @param celltype_col Metadata cell type column.
#' @return List of model results and binned trend data.
#' @export
find_hotspot_distance_genes <- function(seurat_obj,
                                        hotspot_result,
                                        assay = "Xenium",
                                        focus = c("perihotspot", "intrahotspot", "both"),
                                        min_cells_expressing = 10,
                                        n_distance_bins = 10,
                                        max_dist_quantile = 0.95,
                                        p_adj_threshold = 0.05,
                                        exclude_tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
                                        celltype_col = "singleR.predicted.id.brcaAtlas") {
  focus <- match.arg(focus)

  hs_meta <- hotspot_result$tcell_meta
  if (!all(hs_meta$cell_id == colnames(seurat_obj))) {
    stop("Hotspot metadata cell IDs do not match Seurat object cells.")
  }

  seurat_obj$dist_to_hotspot <- hs_meta$dist_to_hotspot
  seurat_obj$signed_dist_hotspot <- hs_meta$signed_dist_hotspot

  if (!is.null(exclude_tcell_labels) && length(exclude_tcell_labels) > 0) {
    cells_keep <- colnames(seurat_obj)[
      !as.character(seurat_obj@meta.data[[celltype_col]]) %in% exclude_tcell_labels
    ]
    seurat_obj <- subset(seurat_obj, cells = cells_keep)
  }

  sdh <- seurat_obj$signed_dist_hotspot
  has_dist <- !is.na(sdh)
  if (focus == "perihotspot") {
    sel <- has_dist & sdh > 0
    cells_use <- colnames(seurat_obj)[sel]
    dist_vec <- sdh[sel]
  } else if (focus == "intrahotspot") {
    sel <- has_dist & sdh <= 0
    cells_use <- colnames(seurat_obj)[sel]
    dist_vec <- abs(sdh[sel])
  } else {
    cells_use <- colnames(seurat_obj)[has_dist]
    dist_vec <- seurat_obj$dist_to_hotspot[has_dist]
  }

  if (length(cells_use) < 20) {
    warning(sprintf(
      "Only %d %s cells (need >= 20); returning NULL. signed_dist_hotspot: %d > 0, %d <= 0, %d NA.",
      length(cells_use), focus, sum(sdh > 0, na.rm = TRUE),
      sum(sdh <= 0, na.rm = TRUE), sum(!has_dist)
    ), call. = FALSE)
    return(NULL)
  }

  sub_obj <- subset(seurat_obj, cells = cells_use)
  dist_vec <- dist_vec[names(dist_vec) %in% colnames(sub_obj)]
  dist_vec <- dist_vec[colnames(sub_obj)]

  max_d <- stats::quantile(dist_vec, max_dist_quantile, na.rm = TRUE)
  keep_idx <- dist_vec <= max_d
  sub_obj <- subset(sub_obj, cells = colnames(sub_obj)[keep_idx])
  dist_vec <- dist_vec[keep_idx]

  expr_mat <- get_data_layer(sub_obj, assay)

  n_expressing <- Matrix::rowSums(expr_mat > 0)
  genes_use <- names(n_expressing[n_expressing >= min_cells_expressing])
  if (length(genes_use) < 2) {
    warn_too_few_genes(n_expressing, min_cells_expressing, ncol(expr_mat), assay)
    return(NULL)
  }
  expr_mat <- expr_mat[genes_use, , drop = FALSE]

  results_lm <- data.frame(
    gene = genes_use,
    coef = NA_real_,
    se = NA_real_,
    t_stat = NA_real_,
    p_value = NA_real_,
    r_squared = NA_real_,
    mean_expr = NA_real_,
    pct_expressing = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes_use)) {
    gene <- genes_use[i]
    y <- as.numeric(expr_mat[gene, ])
    fit <- tryCatch(summary(stats::lm(y ~ dist_vec)), error = function(e) NULL)
    if (!is.null(fit)) {
      results_lm$coef[i] <- fit$coefficients[2, 1]
      results_lm$se[i] <- fit$coefficients[2, 2]
      results_lm$t_stat[i] <- fit$coefficients[2, 3]
      results_lm$p_value[i] <- fit$coefficients[2, 4]
      results_lm$r_squared[i] <- fit$r.squared
    }
    results_lm$mean_expr[i] <- mean(y)
    results_lm$pct_expressing[i] <- mean(y > 0) * 100
  }

  results_lm$p_adj <- stats::p.adjust(results_lm$p_value, method = "BH")
  results_lm$direction <- ifelse(results_lm$coef > 0, "increases_with_distance", "decreases_with_distance")
  results_lm <- results_lm %>% dplyr::arrange(p_adj)

  bin_breaks <- unique(stats::quantile(
    dist_vec,
    probs = seq(0, 1, length.out = n_distance_bins + 1),
    na.rm = TRUE
  ))
  dist_bins <- cut(dist_vec, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)
  bin_midpts <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2

  sig_genes <- results_lm$gene[results_lm$p_adj < p_adj_threshold]
  binned_expr <- NULL
  if (length(sig_genes) > 0) {
    top_genes <- utils::head(sig_genes, min(50, length(sig_genes)))
    binned_expr <- do.call(rbind, lapply(top_genes, function(gene) {
      y <- as.numeric(expr_mat[gene, ])
      bvals <- tapply(y, dist_bins, mean)
      data.frame(
        gene = gene,
        bin = seq_along(bvals),
        bin_midpoint = bin_midpts[seq_along(bvals)],
        mean_expr = as.numeric(bvals),
        stringsAsFactors = FALSE
      )
    }))
  }

  spearman_res <- data.frame(gene = genes_use, rho = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(genes_use)) {
    gene <- genes_use[i]
    y <- as.numeric(expr_mat[gene, ])
    ct <- tryCatch(stats::cor.test(dist_vec, y, method = "spearman", exact = FALSE), error = function(e) NULL)
    if (!is.null(ct)) {
      spearman_res$rho[i] <- unname(ct$estimate)
      spearman_res$p_value[i] <- ct$p.value
    }
  }
  spearman_res$p_adj <- stats::p.adjust(spearman_res$p_value, method = "BH")

  results <- merge(results_lm, spearman_res, by = "gene", suffixes = c("_lm", "_spearman")) %>%
    dplyr::arrange(p_adj_lm)

  list(results = results, binned_expr = binned_expr, n_sig = length(sig_genes), focus = focus)
}

#' Plot T Cell Hotspots on Spatial Coordinates
#'
#' @param seurat_obj Seurat object.
#' @param hotspot_result Output from `define_tcell_hotspots`.
#' @param celltype_col Metadata cell type column.
#' @param tcell_labels T cell labels to highlight.
#' @param sample_label Optional label for title.
#' @param fov_name Field of view/image name for coordinate extraction.
#' @return ggplot object.
#' @export
plot_tcell_hotspots <- function(seurat_obj,
                                hotspot_result,
                                celltype_col = "singleR.predicted.id.brcaAtlas",
                                tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
                                sample_label = "",
                                fov_name = "fov") {
  coords <- get_cell_coords(seurat_obj, fov_name)
  coords$celltype <- as.character(seurat_obj@meta.data[coords$cell_id, celltype_col])
  coords$is_tcell <- coords$celltype %in% tcell_labels

  hs_sf <- hotspot_result$hotspot_sf
  hs_df <- do.call(rbind, lapply(seq_len(nrow(hs_sf)), function(i) {
    poly <- hs_sf$geometry[[i]]
    ring <- as.data.frame(poly[[1]])
    colnames(ring) <- c("x", "y")
    ring$hotspot_id <- i
    ring
  }))

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = coords[!coords$is_tcell, , drop = FALSE],
      ggplot2::aes(x = x, y = y),
      color = "grey80",
      size = 0.3,
      alpha = 0.4
    ) +
    ggplot2::geom_point(
      data = coords[coords$is_tcell, , drop = FALSE],
      ggplot2::aes(x = x, y = y, color = celltype),
      size = 0.8,
      alpha = 0.7
    ) +
    ggplot2::geom_polygon(
      data = hs_df,
      ggplot2::aes(x = x, y = y, group = hotspot_id),
      fill = NA,
      color = "#E74C3C",
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(values = c("T Cells CD4+" = "#3498DB", "T Cells CD8+" = "#E67E22"), name = "T cell type") +
    ggplot2::labs(
      title = ifelse(
        nchar(sample_label) > 0,
        sprintf("T Cell Hotspots - %s", sample_label),
        "T Cell Hotspots"
      ),
      subtitle = sprintf("%d hotspot(s) identified", hotspot_result$n_hotspots),
      x = "x",
      y = "y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(aspect.ratio = 1, legend.position = "right")
}
