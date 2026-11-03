#' Default TLS Signature Genes
#'
#' Returns the default gene set used to score tertiary lymphoid structures
#' (TLS) with [add_tls_score()]. Genes absent from a panel are dropped at
#' scoring time.
#'
#' @return Character vector of gene symbols.
#' @export
#'
#' @examples
#' tls_signature_genes()
tls_signature_genes <- function() {
  c("MS4A1", "CD19", "CD22", "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD74",
    "CD79A", "IL7R", "ITGAE", "CD1D", "CD52", "CD79B", "FCER2", "CR2",
    "PDCD1", "PTGDS", "TRBC2", "CXCL13", "CXCR5")
}

#' Call Tertiary Lymphoid Structures (TLS)
#'
#' Clusters lymphoid (B and T) cells with DBSCAN to find candidate immune
#' aggregates, assigns every cell within `eps` of a clustered lymphoid cell to
#' that aggregate, and calls an aggregate a TLS when it contains enough
#' co-localized B and T cells. Aggregates containing mature dendritic cells
#' (e.g. LAMP3+/DC-LAMP+) are flagged as mature-like.
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param fov_name Field of view/image name (default: "fov").
#' @param celltype_col Metadata column with cell type labels.
#' @param b_types Cell type labels counted as B cells.
#' @param t_types Cell type labels counted as T cells.
#' @param dc_types Cell type labels counted as mature dendritic cells.
#' @param eps DBSCAN neighborhood radius in coordinate units (microns for
#'   Xenium). Also the maximum distance for assigning non-lymphoid cells to an
#'   aggregate.
#' @param minPts DBSCAN minimum lymphoid cells to seed a cluster.
#' @param min_B Minimum B cells for an aggregate to be called a TLS.
#' @param min_T Minimum T cells for an aggregate to be called a TLS.
#' @param min_frac_B Minimum fraction of B cells in the aggregate (excludes
#'   pure T cell infiltrates).
#'
#' @return Seurat object with `TLS_id` (aggregate ID for cells in called TLS,
#'   otherwise `NA`) and `in_TLS` metadata. The per-aggregate summary is stored
#'   in `seurat_obj@misc$tls[[fov_name]]`; retrieve it with [get_tls_summary()].
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- call_tls(seurat_obj, celltype_col = "first_type")
#' get_tls_summary(seurat_obj)
#' }
call_tls <- function(seurat_obj,
                     fov_name = "fov",
                     celltype_col = "first_type",
                     b_types = c("B_Cells"),
                     t_types = c("CD4+_T_Cells", "CD8+_T_Cells"),
                     dc_types = c("LAMP3+_DCs", "IRF7+_DCs"),
                     eps = 50,
                     minPts = 20,
                     min_B = 20,
                     min_T = 20,
                     min_frac_B = 0.05) {
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", celltype_col, "' not found in metadata.", call. = FALSE)
  }

  df <- get_cell_coords(seurat_obj, fov_name)
  df$celltype <- as.character(seurat_obj@meta.data[df$cell_id, celltype_col])
  df$is_B <- df$celltype %in% b_types
  df$is_T <- df$celltype %in% t_types
  df$is_DC <- df$celltype %in% dc_types
  df$is_lymphoid <- df$is_B | df$is_T

  lymph <- df[df$is_lymphoid, , drop = FALSE]
  if (nrow(lymph) == 0) {
    stop("No cells of b_types/t_types found in '", celltype_col, "'.", call. = FALSE)
  }

  db <- dbscan::dbscan(as.matrix(lymph[, c("x", "y")]), eps = eps, minPts = minPts)
  lymph_clustered <- lymph[db$cluster > 0, , drop = FALSE]
  lymph_clustered$cluster <- db$cluster[db$cluster > 0]
  message("Candidate lymphoid aggregates found: ",
          length(unique(lymph_clustered$cluster)))

  df$assigned_cluster <- 0L
  if (nrow(lymph_clustered) > 0) {
    nn <- dbscan::kNN(as.matrix(lymph_clustered[, c("x", "y")]), k = 1,
                      query = as.matrix(df[, c("x", "y")]))
    df$assigned_cluster <- ifelse(nn$dist[, 1] <= eps,
                                  lymph_clustered$cluster[nn$id[, 1]], 0L)
  }

  comp <- df %>%
    dplyr::filter(.data$assigned_cluster > 0) %>%
    dplyr::group_by(cluster = .data$assigned_cluster) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_B = sum(.data$is_B),
      n_T = sum(.data$is_T),
      n_DC = sum(.data$is_DC),
      frac_B = n_B / n_total,
      frac_T = n_T / n_total,
      frac_immune = mean(.data$is_lymphoid | .data$is_DC),
      cx = mean(.data$x), cy = mean(.data$y),
      area_approx = (max(.data$x) - min(.data$x)) * (max(.data$y) - min(.data$y)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      is_TLS = n_B >= min_B & n_T >= min_T & frac_B >= min_frac_B,
      has_matureDC = n_DC > 0,
      maturity = dplyr::case_when(
        is_TLS & has_matureDC ~ "mature-like (DC-LAMP+)",
        is_TLS ~ "early/aggregate",
        TRUE ~ "not TLS"
      )
    )

  tls_ids <- comp$cluster[comp$is_TLS]
  message("TLS called: ", length(tls_ids))

  tls_id <- ifelse(df$assigned_cluster %in% tls_ids, df$assigned_cluster, NA_integer_)
  seurat_obj$TLS_id <- expand_to_all_cells(tls_id, df$cell_id, seurat_obj, "TLS_id")
  seurat_obj$in_TLS <- expand_to_all_cells(!is.na(tls_id), df$cell_id, seurat_obj, "in_TLS")

  if (is.null(seurat_obj@misc$tls)) seurat_obj@misc$tls <- list()
  seurat_obj@misc$tls[[fov_name]] <- as.data.frame(comp)
  seurat_obj
}

#' Get TLS Aggregate Summary
#'
#' Retrieves the per-aggregate composition table stored by [call_tls()].
#'
#' @param seurat_obj Seurat object processed with [call_tls()].
#' @param fov_name Field of view/image name (default: "fov").
#' @param called_only If `TRUE` (default), return only aggregates called as TLS.
#' @return Data frame with one row per aggregate: cell counts (`n_total`,
#'   `n_B`, `n_T`, `n_DC`), fractions, centroid (`cx`, `cy`), `is_TLS` and
#'   `maturity`. `NULL` with a warning if [call_tls()] has not been run.
#' @export
get_tls_summary <- function(seurat_obj, fov_name = "fov", called_only = TRUE) {
  comp <- seurat_obj@misc$tls[[fov_name]]
  if (is.null(comp)) {
    warning("No TLS calls stored for FOV '", fov_name, "'. Run call_tls() first.",
            call. = FALSE)
    return(NULL)
  }
  if (called_only) comp <- comp[comp$is_TLS, , drop = FALSE]
  rownames(comp) <- NULL
  comp
}

#' Add a TLS Signature Score
#'
#' Scores each cell for a TLS gene signature with [Seurat::AddModuleScore()],
#' which corrects for background expression using control gene sets. Genes not
#' in the assay are dropped with a message.
#'
#' @param seurat_obj Seurat object.
#' @param genes Signature genes (default: [tls_signature_genes()]).
#' @param assay Assay name. If its `data` layer is empty,
#'   [Seurat::NormalizeData()] is run with default settings.
#' @param score_name Metadata column to store the score in.
#' @param ctrl Number of control genes per signature gene. Defaults to the
#'   number of genes present, bounded to 5-50 and capped by the genes
#'   available per expression bin, which suits small panels.
#' @param nbin Number of expression bins for control gene selection; lower
#'   this if AddModuleScore fails on small panels.
#' @param min_genes Minimum signature genes that must be present.
#' @param seed Random seed passed to AddModuleScore.
#' @return Seurat object with the score in metadata column `score_name`.
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- add_tls_score(seurat_obj, assay = "Xenium")
#' }
add_tls_score <- function(seurat_obj,
                          genes = tls_signature_genes(),
                          assay = "Xenium",
                          score_name = "TLS_score",
                          ctrl = NULL,
                          nbin = 24,
                          min_genes = 3,
                          seed = 1) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found. Available assays: ",
         paste(names(seurat_obj@assays), collapse = ", "), call. = FALSE)
  }
  present <- intersect(genes, rownames(seurat_obj[[assay]]))
  missing <- setdiff(genes, present)
  if (length(missing) > 0) {
    message("Dropped ", length(missing), " gene(s) not in assay '", assay, "': ",
            paste(missing, collapse = ", "))
  }
  if (length(present) < min_genes) {
    stop("Only ", length(present), " signature gene(s) present in assay '", assay,
         "' (need >= ", min_genes, ").", call. = FALSE)
  }
  message("Scoring with ", length(present), " gene(s): ", paste(present, collapse = ", "))

  data_mat <- suppressWarnings(tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data"),
    error = function(e) NULL
  ))
  if (is.null(data_mat) || nrow(data_mat) == 0 || ncol(data_mat) == 0) {
    message("The 'data' layer of assay '", assay, "' is empty; running ",
            "Seurat::NormalizeData() with default settings.")
    seurat_obj <- Seurat::NormalizeData(seurat_obj, assay = assay, verbose = FALSE)
  }

  if (is.null(ctrl)) {
    # Control genes are sampled without replacement from each expression bin,
    # so keep ctrl below the genes available per bin on small panels.
    per_bin <- max(1, floor(nrow(seurat_obj[[assay]]) / nbin) - 1)
    ctrl <- min(50, max(5, length(present)), per_bin)
  }
  tmp_name <- paste0(score_name, "_tmp")
  seurat_obj <- tryCatch(
    Seurat::AddModuleScore(
      seurat_obj, features = list(present), assay = assay, name = tmp_name,
      ctrl = ctrl, nbin = nbin, seed = seed
    ),
    error = function(e) {
      stop("AddModuleScore failed (", conditionMessage(e), "). Try lowering ",
           "`nbin` and/or `ctrl`.", call. = FALSE)
    }
  )
  seurat_obj@meta.data[[score_name]] <- seurat_obj@meta.data[[paste0(tmp_name, "1")]]
  seurat_obj@meta.data[[paste0(tmp_name, "1")]] <- NULL
  seurat_obj
}

#' Summarize TLS Signature Scores
#'
#' Sanity check that called TLS score higher than the rest of the tissue.
#'
#' @param seurat_obj Seurat object with [call_tls()] and [add_tls_score()]
#'   results.
#' @param score_col Metadata column with the signature score.
#' @return List with `by_status` (in vs. outside TLS), `by_tls` (per TLS,
#'   ordered by mean score) and `wilcox_p` (two-sided Wilcoxon test p-value).
#' @export
summarize_tls_score <- function(seurat_obj, score_col = "TLS_score") {
  df <- tls_score_df(seurat_obj, score_col)
  df <- df[!is.na(df$in_TLS), , drop = FALSE]

  summ <- function(d) dplyr::summarise(
    d,
    mean_score = mean(.data$score, na.rm = TRUE),
    median_score = stats::median(.data$score, na.rm = TRUE),
    n = dplyr::n(), .groups = "drop"
  )

  by_status <- summ(dplyr::group_by(df, in_TLS))
  by_tls <- df %>%
    dplyr::filter(!is.na(.data$TLS_id)) %>%
    dplyr::group_by(TLS_id) %>%
    summ() %>%
    dplyr::arrange(.data$mean_score)

  list(by_status = as.data.frame(by_status),
       by_tls = as.data.frame(by_tls),
       wilcox_p = tls_wilcox_p(df))
}

#' Plot Called TLS on the Tissue
#'
#' Spatial map of all cells (grey) with cells in called TLS colored by TLS ID,
#' B cells within TLS in red, and TLS IDs labeled at their centroids.
#'
#' @param seurat_obj Seurat object processed with [call_tls()].
#' @param fov_name Field of view/image name (default: "fov").
#' @param celltype_col Metadata column with cell type labels.
#' @param b_types Cell type labels to highlight as B cells.
#' @param point_size Point size for TLS cells.
#' @return ggplot object.
#' @export
plot_tls_map <- function(seurat_obj,
                         fov_name = "fov",
                         celltype_col = "first_type",
                         b_types = c("B_Cells"),
                         point_size = 0.5) {
  df <- tls_cell_df(seurat_obj, fov_name)
  df$is_B <- as.character(seurat_obj@meta.data[df$cell_id, celltype_col]) %in% b_types
  tls_df <- df[df$in_TLS %in% TRUE, , drop = FALSE]
  centroids <- tls_centroids(tls_df)

  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(color = "grey85", size = 0.2) +
    ggplot2::geom_point(data = tls_df, ggplot2::aes(color = factor(.data$TLS_id)),
                        size = point_size) +
    ggplot2::geom_point(data = tls_df[tls_df$is_B, , drop = FALSE],
                        color = "red", size = point_size) +
    ggplot2::geom_text(data = centroids,
                       ggplot2::aes(.data$cx, .data$cy, label = .data$TLS_id),
                       size = 3, fontface = "bold") +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::labs(title = "Called TLS (B cells in red)", color = "TLS id")
}

#' Plot TLS Cellular Composition
#'
#' Stacked bar chart of cell type fractions within each called TLS.
#'
#' @param seurat_obj Seurat object processed with [call_tls()].
#' @param celltype_col Metadata column with cell type labels.
#' @return ggplot object.
#' @export
plot_tls_composition <- function(seurat_obj, celltype_col = "first_type") {
  meta <- seurat_obj@meta.data
  if (!"TLS_id" %in% colnames(meta)) {
    stop("Metadata column 'TLS_id' not found. Run call_tls() first.", call. = FALSE)
  }
  df <- data.frame(TLS_id = meta$TLS_id,
                   celltype = as.character(meta[[celltype_col]]),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$TLS_id), , drop = FALSE]
  if (nrow(df) == 0) stop("No cells in called TLS.", call. = FALSE)

  counts <- dplyr::count(df, .data$TLS_id, .data$celltype)
  ggplot2::ggplot(counts, ggplot2::aes(factor(.data$TLS_id), .data$n, fill = .data$celltype)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::scale_fill_manual(values = celltype_palette(unique(counts$celltype))) +
    ggplot2::labs(x = "TLS id", y = "fraction", fill = "Cell type",
                  title = "TLS cellular composition") +
    ggplot2::theme_minimal()
}

#' Plot TLS Signature Score Over Space
#'
#' Colors all cells by TLS signature score (clipped to the 1st-99th percentile)
#' on a dark tissue background, optionally highlighting called TLS.
#'
#' @param seurat_obj Seurat object with [add_tls_score()] results (and
#'   [call_tls()] results unless `highlight = "none"`).
#' @param fov_name Field of view/image name (default: "fov").
#' @param score_col Metadata column with the signature score.
#' @param highlight How to mark called TLS: `"outline"` draws a white concave
#'   hull perimeter around each TLS, `"ring"` circles each TLS cell in cyan,
#'   `"none"` shows the score only.
#' @param label_tls Label each TLS with its ID at its centroid.
#' @param concavity Concave hull concavity for `"outline"` (lower = tighter).
#' @param point_size Point size for cells.
#' @return ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_tls_score(seurat_obj, highlight = "outline")
#' }
plot_tls_score <- function(seurat_obj,
                           fov_name = "fov",
                           score_col = "TLS_score",
                           highlight = c("outline", "ring", "none"),
                           label_tls = FALSE,
                           concavity = 2,
                           point_size = 0.25) {
  highlight <- match.arg(highlight)
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column '", score_col, "' not found. Run add_tls_score() first.",
         call. = FALSE)
  }
  df <- tls_cell_df(seurat_obj, fov_name, require_tls = highlight != "none" || label_tls)
  df$score <- seurat_obj@meta.data[df$cell_id, score_col]
  q <- stats::quantile(df$score, c(0.01, 0.99), na.rm = TRUE)
  df$score_clip <- pmin(pmax(df$score, q[1]), q[2])

  tissue_hull <- df[grDevices::chull(df$x, df$y), c("x", "y")]
  tls_df <- df[df$in_TLS %in% TRUE, , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_polygon(data = tissue_hull, fill = "grey15", color = NA) +
    ggplot2::geom_point(ggplot2::aes(color = .data$score_clip), size = point_size) +
    ggplot2::scale_color_viridis_c(option = "magma", name = "TLS\nscore")

  if (highlight == "ring" && nrow(tls_df) > 0) {
    p <- p + ggplot2::geom_point(data = tls_df, shape = 21, color = "cyan", fill = NA,
                                 stroke = 0.25, size = 0.1)
  } else if (highlight == "outline" && nrow(tls_df) > 0) {
    hulls <- tls_hulls(tls_df, concavity)
    p <- p + ggplot2::geom_polygon(data = hulls,
                                   ggplot2::aes(.data$x, .data$y, group = .data$TLS_id),
                                   fill = NA, color = "white", linewidth = 0.6)
  }
  if (label_tls && nrow(tls_df) > 0) {
    p <- p + ggplot2::geom_text(data = tls_centroids(tls_df),
                                ggplot2::aes(.data$cx, .data$cy, label = .data$TLS_id),
                                color = "white", fontface = "bold", size = 3)
  }

  title <- switch(highlight,
                  outline = "TLS signature score with perimeter outlines",
                  ring = "TLS signature score with called TLS highlighted (cyan rings)",
                  none = "TLS signature score")
  p + ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = title)
}

#' Compare TLS Signature Score In vs. Outside TLS
#'
#' Plots the per-cell TLS signature score for cells inside called TLS versus
#' the rest of the tissue, annotated with a two-sided Wilcoxon rank-sum test
#' p-value.
#'
#' @param seurat_obj Seurat object with [call_tls()] and [add_tls_score()]
#'   results.
#' @param score_col Metadata column with the signature score.
#' @param type `"beeswarm"` (quasirandom points with boxplot; uses
#'   \pkg{ggbeeswarm} if installed, otherwise jittered points) or `"violin"`
#'   (violin with quartile lines).
#' @return ggplot object.
#' @export
plot_tls_score_comparison <- function(seurat_obj,
                                      score_col = "TLS_score",
                                      type = c("beeswarm", "violin")) {
  type <- match.arg(type)
  df <- tls_score_df(seurat_obj, score_col)
  df <- df[!is.na(df$in_TLS), , drop = FALSE]
  df$group <- factor(ifelse(df$in_TLS, "In TLS", "Outside TLS"),
                     levels = c("Outside TLS", "In TLS"))
  cols <- c("Outside TLS" = "grey50", "In TLS" = "#B63679")

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$group, .data$score))
  if (type == "violin") {
    p <- p +
      ggplot2::geom_violin(ggplot2::aes(fill = .data$group), width = 0.6, trim = TRUE,
                           draw_quantiles = c(0.25, 0.5, 0.75)) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::theme_bw(base_size = 13)
  } else {
    pts <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
      ggbeeswarm::geom_quasirandom(ggplot2::aes(color = .data$group),
                                   size = 0.15, alpha = 0.15, width = 0.4)
    } else {
      ggplot2::geom_jitter(ggplot2::aes(color = .data$group),
                           size = 0.15, alpha = 0.15, width = 0.3, height = 0)
    }
    p <- p + pts +
      ggplot2::geom_boxplot(fill = NA, color = "black", width = 0.25,
                            outlier.shape = NA, linewidth = 0.5) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::theme_minimal(base_size = 13)
  }

  pval <- tls_wilcox_p(df)
  if (!is.na(pval)) {
    y_max <- max(df$score, na.rm = TRUE)
    y_rng <- diff(range(df$score, na.rm = TRUE))
    y_br <- y_max + 0.05 * y_rng
    p_lab <- format.pval(pval, digits = 2)
    p_lab <- if (startsWith(p_lab, "<")) paste("p", p_lab) else paste("p =", p_lab)
    p <- p +
      ggplot2::annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br) +
      ggplot2::annotate("text", x = 1.5, y = y_br + 0.04 * y_rng, label = p_lab)
  }

  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(x = NULL, y = "TLS signature score",
                  title = "TLS signature score: called TLS vs. rest of tissue") +
    ggplot2::theme(legend.position = "none")
}

# Coordinates plus TLS calls for cells in the FOV.
tls_cell_df <- function(seurat_obj, fov_name, require_tls = TRUE) {
  df <- get_cell_coords(seurat_obj, fov_name)
  meta <- seurat_obj@meta.data
  has_tls <- all(c("TLS_id", "in_TLS") %in% colnames(meta))
  if (require_tls && !has_tls) {
    stop("Metadata columns 'TLS_id'/'in_TLS' not found. Run call_tls() first.",
         call. = FALSE)
  }
  df$TLS_id <- if (has_tls) meta[df$cell_id, "TLS_id"] else NA
  df$in_TLS <- if (has_tls) meta[df$cell_id, "in_TLS"] else NA
  df
}

# Score plus TLS calls for all cells.
tls_score_df <- function(seurat_obj, score_col) {
  meta <- seurat_obj@meta.data
  needed <- c(score_col, "TLS_id", "in_TLS")
  absent <- setdiff(needed, colnames(meta))
  if (length(absent) > 0) {
    stop("Metadata column(s) not found: ", paste(absent, collapse = ", "),
         ". Run call_tls() and add_tls_score() first.", call. = FALSE)
  }
  data.frame(score = meta[[score_col]], TLS_id = meta$TLS_id, in_TLS = meta$in_TLS)
}

tls_wilcox_p <- function(df) {
  in_s <- df$score[df$in_TLS %in% TRUE]
  out_s <- df$score[df$in_TLS %in% FALSE]
  if (length(in_s) == 0 || length(out_s) == 0) return(NA_real_)
  suppressWarnings(stats::wilcox.test(in_s, out_s)$p.value)
}

tls_centroids <- function(tls_df) {
  tls_df %>%
    dplyr::group_by(TLS_id) %>%
    dplyr::summarise(cx = mean(.data$x), cy = mean(.data$y), .groups = "drop")
}

# Concave hull perimeter per TLS (TLS with < 3 cells are skipped).
tls_hulls <- function(tls_df, concavity = 2) {
  parts <- lapply(split(tls_df, tls_df$TLS_id), function(d) {
    if (nrow(d) < 3) return(NULL)
    h <- concaveman::concaveman(as.matrix(d[, c("x", "y")]), concavity = concavity)
    data.frame(TLS_id = d$TLS_id[1], x = h[, 1], y = h[, 2])
  })
  do.call(rbind, parts)
}
