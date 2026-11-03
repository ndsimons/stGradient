#' Compute Signed Distance From Tumor Boundary
#'
#' Computes per-cell signed distance to the nearest cancer region polygon, using
#' the same per-region polygons and definition as [calculate_boundary_distances()].
#' `dist_to_boundary` is therefore identical to `boundary_distance` computed from
#' the same polygons, and stroma between separate tumor nests is classed as
#' outside the tumor.
#'
#' @param seurat_obj Seurat object.
#' @param polygons sf polygon object. Defaults to the polygons stored by
#'   [add_dbscan_polygons()], retrieved with [get_cancer_regions()].
#' @param fov_name Field of view/image name for coordinate extraction.
#' @param polygon_col,intratumoral_label,concavity Deprecated and ignored. The
#'   boundary used to be a single concave hull around all intratumoral cells,
#'   which merged separate tumor nests and labelled the stroma between them as
#'   intratumoral.
#' @return Seurat object with `dist_to_boundary` and `abs_dist_to_boundary` metadata.
#' @export
compute_distance_from_boundary <- function(seurat_obj,
                                           polygons = NULL,
                                           fov_name = "fov",
                                           polygon_col = NULL,
                                           intratumoral_label = NULL,
                                           concavity = NULL) {
  deprecated <- c("polygon_col", "intratumoral_label", "concavity")
  used <- deprecated[!vapply(list(polygon_col, intratumoral_label, concavity),
                             is.null, logical(1))]
  if (length(used) > 0) {
    warning("`", paste(used, collapse = "`, `"), "` ",
            if (length(used) == 1) "is" else "are",
            " deprecated and ignored; the boundary is now taken from the cancer ",
            "region polygons (see `polygons`).", call. = FALSE)
  }

  if (is.null(polygons)) {
    polygons <- get_cancer_regions(seurat_obj, fov_name)
  }
  signed_distance <- calculate_boundary_distances(seurat_obj, polygons, fov_name)

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
#' @param max_distance Maximum distance from the boundary to include (for
#'   `intratumoral`, maximum absolute distance). Cells farther away are excluded
#'   from model fitting and binning. Combined with `max_dist_quantile`: the
#'   smaller of the two cutoffs is used. Set to `NULL` to rely on the quantile
#'   alone.
#' @param p_adj_threshold BH adjusted p-value threshold for significant genes.
#' @return List containing model results, binned expression, and metadata.
#' @export
find_distance_associated_genes <- function(seurat_obj,
                                           assay = "Xenium",
                                           min_cells_expressing = 10,
                                           n_distance_bins = 10,
                                           focus = c("peritumoral", "intratumoral", "both"),
                                           max_dist_quantile = 0.95,
                                           max_distance = 200,
                                           p_adj_threshold = 0.05) {
  focus <- match.arg(focus)

  if (!is.null(max_distance) &&
      !(is.numeric(max_distance) && length(max_distance) == 1 &&
        is.finite(max_distance) && max_distance > 0)) {
    stop("`max_distance` must be a single positive number or NULL.", call. = FALSE)
  }

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

  max_dist <- unname(stats::quantile(dist_values, max_dist_quantile, na.rm = TRUE))
  if (!is.null(max_distance)) {
    max_dist <- min(max_dist, max_distance)
  }
  keep <- dist_values <= max_dist
  if (sum(keep) < 20) {
    warning(sprintf(
      "Only %d %s cells within distance %g (need >= 20); returning NULL. Increase max_distance.",
      sum(keep), focus, max_dist
    ), call. = FALSE)
    return(NULL)
  }
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
    n_distance_bins = n_distance_bins,
    max_distance = max_distance
  )
}

#' Plot Distance-Associated Gene Results
#'
#' @param result_list Output list from `find_distance_associated_genes`.
#' @param top_n Maximum number of significant genes shown in trend facets.
#' @param sample_label Optional label for plot titles.
#' @return List with `volcano` and `trends` ggplot objects. The trends x-axis
#'   extends to the `max_distance` used in `find_distance_associated_genes()`.
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

      max_distance <- result_list$max_distance
      if (!is.null(max_distance)) {
        x_min <- min(0, plot_data$bin_midpoint, na.rm = TRUE)
        p_trends <- p_trends + ggplot2::coord_cartesian(xlim = c(x_min, max_distance))
      }
    }
  }

  list(volcano = p_volcano, trends = p_trends)
}
