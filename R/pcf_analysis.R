#' Build an Outside-Polygon Point Pattern for PCF
#'
#' Constructs a marked point pattern (`ppp`) using only cells outside the supplied
#' tumor polygons. This helper is used internally by PCF analysis functions.
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param polygons sf polygon object representing ROIs.
#' @param celltype_col Metadata column containing cell type labels.
#' @param fov_name Field of view/image name passed to `Seurat::GetTissueCoordinates`.
#' @param min_cells Minimum number of cells required per cell type.
#' @param window_obj Optional `spatstat.geom::owin` window. If `NULL`, a bounding-box
#'   window is built from all coordinates.
#' @return A list with `ppp_obj`, `cells_outside`, and `window`.
#' @keywords internal
build_outside_ppp <- function(seurat_obj,
                              polygons,
                              celltype_col = "singleR.predicted.id.brcaAtlas",
                              fov_name = "fov",
                              min_cells = 10,
                              window_obj = NULL) {
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
    stop("Package 'spatstat.geom' is required. Install with: install.packages('spatstat.geom')")
  }

  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", celltype_col, "' not found in seurat_obj@meta.data")
  }

  coords_df <- get_cell_coords(seurat_obj, fov_name)
  check_polygons(polygons)
  coords_df$cellType <- as.character(seurat_obj@meta.data[coords_df$cell_id, celltype_col])

  coords_sf <- sf::st_as_sf(coords_df, coords = c("x", "y"))
  inside_mat <- sf::st_intersects(coords_sf, polygons, sparse = FALSE)
  coords_df$inside_polygon <- apply(inside_mat, 1, any)

  cells_outside <- coords_df[!coords_df$inside_polygon, , drop = FALSE]
  cells_outside <- cells_outside[!is.na(cells_outside$cellType), , drop = FALSE]

  if (nrow(cells_outside) < min_cells) {
    stop("Not enough cells outside polygons to run PCF analysis.")
  }

  cell_type_counts <- table(cells_outside$cellType)
  valid_types <- names(cell_type_counts[cell_type_counts >= min_cells])
  cells_outside <- cells_outside[cells_outside$cellType %in% valid_types, , drop = FALSE]
  cells_outside$cellType <- factor(cells_outside$cellType)

  if (nrow(cells_outside) < min_cells || nlevels(cells_outside$cellType) < 2) {
    stop("Insufficient outside-polygon cells/cell types after min_cells filtering.")
  }

  if (is.null(window_obj)) {
    window_obj <- spatstat.geom::owin(
      xrange = range(coords_df$x, na.rm = TRUE),
      yrange = range(coords_df$y, na.rm = TRUE)
    )
  }

  ppp_obj <- spatstat.geom::ppp(
    x = cells_outside$x,
    y = cells_outside$y,
    marks = cells_outside$cellType,
    window = window_obj,
    checkdup = FALSE
  )

  list(ppp_obj = ppp_obj, cells_outside = cells_outside, window = window_obj)
}

#' Calculate Outside-Polygon Pair Correlation Summary
#'
#' Computes pair correlation function (PCF) values for all requested cell-type pairs
#' using only cells outside the tumor polygons.
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param polygons sf polygon object representing ROIs.
#' @param distance Target distance at which to extract PCF value.
#' @param min_cells Minimum number of cells required per cell type.
#' @param celltype_col Metadata column containing cell type labels.
#' @param fov_name Field of view/image name passed to `Seurat::GetTissueCoordinates`.
#' @param include_self_pairs Logical; include same-type pairs.
#' @param ordered_pairs Logical; if `TRUE`, computes ordered pairs (A,B and B,A).
#'   If `FALSE`, computes each pair once.
#' @param window_obj Optional `spatstat.geom::owin` window.
#' @param show_progress Logical; show a text progress bar.
#' @return Data frame with columns `type1`, `type2`, `g_value`, `n_type1`, `n_type2`.
#' @export
calculate_pcf_outside_summary <- function(seurat_obj,
                                          polygons,
                                          distance = 30,
                                          min_cells = 10,
                                          celltype_col = "singleR.predicted.id.brcaAtlas",
                                          fov_name = "fov",
                                          include_self_pairs = TRUE,
                                          ordered_pairs = FALSE,
                                          window_obj = NULL,
                                          show_progress = interactive()) {
  if (!requireNamespace("spatstat.explore", quietly = TRUE)) {
    stop("Package 'spatstat.explore' is required. Install with: install.packages('spatstat.explore')")
  }

  prep <- build_outside_ppp(
    seurat_obj = seurat_obj,
    polygons = polygons,
    celltype_col = celltype_col,
    fov_name = fov_name,
    min_cells = min_cells,
    window_obj = window_obj
  )

  ppp_obj <- prep$ppp_obj
  marks_vec <- as.character(spatstat.geom::marks(ppp_obj))
  cell_types <- levels(spatstat.geom::marks(ppp_obj))

  if (ordered_pairs) {
    pair_grid <- expand.grid(type1 = cell_types, type2 = cell_types, stringsAsFactors = FALSE)
    if (!include_self_pairs) {
      pair_grid <- pair_grid[pair_grid$type1 != pair_grid$type2, , drop = FALSE]
    }
  } else {
    pair_list <- utils::combn(cell_types, 2, simplify = FALSE)
    pair_grid <- do.call(rbind, lapply(pair_list, function(x) {
      data.frame(type1 = x[1], type2 = x[2], stringsAsFactors = FALSE)
    }))
    if (include_self_pairs) {
      self_df <- data.frame(type1 = cell_types, type2 = cell_types, stringsAsFactors = FALSE)
      pair_grid <- rbind(self_df, pair_grid)
    }
    if (nrow(pair_grid) == 0) {
      pair_grid <- data.frame(type1 = character(0), type2 = character(0), stringsAsFactors = FALSE)
    }
  }

  if (nrow(pair_grid) == 0) {
    return(data.frame(
      type1 = character(0),
      type2 = character(0),
      g_value = numeric(0),
      n_type1 = integer(0),
      n_type2 = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  n_type <- table(marks_vec)
  out_rows <- vector("list", nrow(pair_grid))

  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(pair_grid), style = 3)
    on.exit({
      close(pb)
      cat("\n")
    }, add = TRUE)
  }

  for (k in seq_len(nrow(pair_grid))) {
    type1 <- pair_grid$type1[k]
    type2 <- pair_grid$type2[k]

    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, k)
    }

    out_rows[[k]] <- tryCatch({
      pcf_ij <- spatstat.explore::pcfcross(ppp_obj, i = type1, j = type2)

      idx <- which.min(abs(pcf_ij$r - distance))
      g_val <- if (length(idx) > 0) pcf_ij$iso[idx] else NA_real_

      data.frame(
        type1 = type1,
        type2 = type2,
        g_value = g_val,
        n_type1 = as.integer(n_type[type1]),
        n_type2 = as.integer(n_type[type2]),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        type1 = type1,
        type2 = type2,
        g_value = NA_real_,
        n_type1 = as.integer(n_type[type1]),
        n_type2 = as.integer(n_type[type2]),
        stringsAsFactors = FALSE
      )
    })
  }

  do.call(rbind, out_rows)
}

#' Plot Outside-Polygon PCF Heatmap
#'
#' Converts PCF summary output into a symmetric matrix and visualizes it as a heatmap.
#'
#' @param pcf_summary Output from `calculate_pcf_outside_summary`.
#' @param midpoint Midpoint value for diverging color scale.
#' @return ggplot object.
#' @export
plot_pcf_outside_heatmap <- function(pcf_summary, midpoint = 1) {
  if (nrow(pcf_summary) == 0) {
    stop("pcf_summary is empty.")
  }

  cell_types <- sort(unique(c(pcf_summary$type1, pcf_summary$type2)))
  pcf_mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(cell_types))
  rownames(pcf_mat) <- cell_types
  colnames(pcf_mat) <- cell_types

  for (i in seq_len(nrow(pcf_summary))) {
    t1 <- pcf_summary$type1[i]
    t2 <- pcf_summary$type2[i]
    v <- pcf_summary$g_value[i]
    pcf_mat[t1, t2] <- v
    pcf_mat[t2, t1] <- v
  }

  pcf_df <- as.data.frame(as.table(pcf_mat), stringsAsFactors = FALSE)
  colnames(pcf_df) <- c("type1", "type2", "g_value")

  .data <- NULL

  ggplot2::ggplot(
    pcf_df,
    ggplot2::aes(x = .data$type2, y = .data$type1, fill = .data$g_value)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "steelblue",
      mid = "white",
      high = "firebrick",
      midpoint = midpoint,
      na.value = "grey90"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Outside-Polygon Pair Correlation Heatmap",
      x = "Cell type",
      y = "Cell type",
      fill = "g(r)"
    )
}

#' Plot Outside-Polygon PCF Envelopes for All Cell-Type Pairs
#'
#' Computes and plots envelope-based pair correlation curves for all requested cell-type
#' pairs outside polygons. When `file` is provided, each pair is saved as its own PNG
#' named `<file stem>_<type1>_vs_<type2>.png`.
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param polygons sf polygon object representing ROIs.
#' @param celltype_col Metadata column containing cell type labels.
#' @param fov_name Field of view/image name passed to `Seurat::GetTissueCoordinates`.
#' @param min_cells Minimum number of cells required per cell type.
#' @param include_self_pairs Logical; include same-type pairs.
#' @param ordered_pairs Logical; if `TRUE`, computes ordered pairs.
#' @param nsim Number of simulations for envelopes.
#' @param bw Bandwidth passed to `pcfcross` in envelope evaluation.
#' @param window_obj Optional `spatstat.geom::owin` window.
#' @param file Optional base PNG file path (e.g. `"pcf_envelopes.png"`). One PNG is
#'   written per pair, with the pair name appended to the file stem.
#' @param width,height PNG size in inches.
#' @param res PNG resolution in pixels per inch.
#' @param verbose Logical; print progress messages.
#' @return Named list of envelope objects.
#' @export
plot_pcf_outside_pair_envelopes <- function(seurat_obj,
                                            polygons,
                                            celltype_col = "singleR.predicted.id.brcaAtlas",
                                            fov_name = "fov",
                                            min_cells = 10,
                                            include_self_pairs = FALSE,
                                            ordered_pairs = FALSE,
                                            nsim = 10,
                                            bw = 10,
                                            window_obj = NULL,
                                            file = NULL,
                                            width = 7,
                                            height = 5,
                                            res = 150,
                                            verbose = TRUE) {
  if (!requireNamespace("spatstat.explore", quietly = TRUE)) {
    stop("Package 'spatstat.explore' is required. Install with: install.packages('spatstat.explore')")
  }

  prep <- build_outside_ppp(
    seurat_obj = seurat_obj,
    polygons = polygons,
    celltype_col = celltype_col,
    fov_name = fov_name,
    min_cells = min_cells,
    window_obj = window_obj
  )

  ppp_obj <- prep$ppp_obj
  cell_types <- levels(spatstat.geom::marks(ppp_obj))

  if (ordered_pairs) {
    pair_grid <- expand.grid(type1 = cell_types, type2 = cell_types, stringsAsFactors = FALSE)
    if (!include_self_pairs) {
      pair_grid <- pair_grid[pair_grid$type1 != pair_grid$type2, , drop = FALSE]
    }
  } else {
    pair_list <- utils::combn(cell_types, 2, simplify = FALSE)
    pair_grid <- do.call(rbind, lapply(pair_list, function(x) {
      data.frame(type1 = x[1], type2 = x[2], stringsAsFactors = FALSE)
    }))
    if (include_self_pairs) {
      pair_grid <- rbind(
        data.frame(type1 = cell_types, type2 = cell_types, stringsAsFactors = FALSE),
        pair_grid
      )
    }
  }

  if (nrow(pair_grid) == 0) {
    return(list())
  }

  if (!is.null(file)) {
    file_stem <- sub("\\.png$", "", file, ignore.case = TRUE)
    pair_files <- paste0(
      file_stem, "_",
      gsub("[^A-Za-z0-9.+-]+", "_", paste(pair_grid$type1, "vs", pair_grid$type2)),
      ".png"
    )
  }

  env_list <- vector("list", nrow(pair_grid))
  for (k in seq_len(nrow(pair_grid))) {
    type1 <- pair_grid$type1[k]
    type2 <- pair_grid$type2[k]
    pair_name <- paste(type1, "vs", type2)

    if (isTRUE(verbose)) {
      message("Computing envelope: ", pair_name, " (", k, "/", nrow(pair_grid), ")")
    }

    png_open <- FALSE
    env_list[[k]] <- tryCatch({
      env_obj <- spatstat.explore::envelope(
        ppp_obj,
        fun = spatstat.explore::pcfcross,
        i = type1,
        j = type2,
        nsim = nsim,
        bw = bw,
        verbose = FALSE
      )

      if (!is.null(file)) {
        grDevices::png(pair_files[k], width = width, height = height, units = "in", res = res)
        png_open <- TRUE
      }

      plot(
        env_obj,
        main = pair_name,
        xlab = "Distance r (microns)",
        ylab = "g(r)",
        shade = c("hi", "lo"),
        legend = FALSE
      )

      env_obj
    }, error = function(e) {
      if (isTRUE(verbose)) {
        message("Skipping ", pair_name, ": ", e$message)
      }
      NULL
    }, finally = if (png_open) grDevices::dev.off())
  }

  names(env_list) <- paste(pair_grid$type1, "vs", pair_grid$type2)
  env_list
}

#' Plot a Single Outside-Polygon PCF Pair as a ggplot
#'
#' Plots the pair correlation function (with simulation envelopes) for one
#' cell-type pair from the list returned by [plot_pcf_outside_pair_envelopes()],
#' as a ggplot object suitable for embedding in R Markdown documents.
#'
#' @param env_list Named list of envelope objects returned by
#'   [plot_pcf_outside_pair_envelopes()].
#' @param celltype1 First cell type label (character).
#' @param celltype2 Second cell type label (character). If
#'   `"<celltype1> vs <celltype2>"` is not in `env_list`, the reversed pair is
#'   used, since unordered envelope lists contain only one direction.
#' @return A ggplot object with the observed g(r) curve and simulation envelope bands.
#' @export
#'
#' @examples
#' \dontrun{
#' env_list <- plot_pcf_outside_pair_envelopes(seurat_obj, auto_polygons)
#' plot_pcf_pair_gg(env_list, "TypeA", "TypeB")
#' }
plot_pcf_pair_gg <- function(env_list, celltype1, celltype2) {
  if (!is.list(env_list) || inherits(env_list, "envelope") || is.null(names(env_list))) {
    stop("`env_list` must be the named list returned by plot_pcf_outside_pair_envelopes().",
         call. = FALSE)
  }

  pair_name <- paste(celltype1, "vs", celltype2)
  reversed_name <- paste(celltype2, "vs", celltype1)
  if (pair_name %in% names(env_list)) {
    pair <- c(celltype1, celltype2)
  } else if (reversed_name %in% names(env_list)) {
    pair_name <- reversed_name
    pair <- c(celltype2, celltype1)
  } else {
    stop("Pair '", celltype1, " vs ", celltype2, "' not found in env_list. Available pairs:\n  ",
         paste(names(env_list), collapse = "\n  "), call. = FALSE)
  }

  env_obj <- env_list[[pair_name]]
  if (is.null(env_obj)) {
    stop("Envelope for '", pair_name, "' is NULL; it failed during ",
         "plot_pcf_outside_pair_envelopes() (rerun with verbose = TRUE to see why).",
         call. = FALSE)
  }

  nsim <- attr(env_obj, "einfo")$nsim
  subtitle <- if (is.null(nsim)) {
    "Outside-polygon PCF"
  } else {
    sprintf("Outside-polygon PCF  |  envelope from %d simulations", as.integer(nsim))
  }
  env_df <- as.data.frame(env_obj)

  ggplot2::ggplot(env_df, ggplot2::aes(x = .data$r)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      fill = "steelblue", alpha = 0.25
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(ggplot2::aes(y = .data$obs), color = "steelblue", linewidth = 1) +
    ggplot2::labs(
      title = paste0(pair[1], "  \u2194  ", pair[2]),
      subtitle = subtitle,
      x = "Distance r (\u03bcm)",
      y = "g(r)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
