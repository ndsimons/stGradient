#' Build the Outside-Tumor Observation Window for PCF
#'
#' Builds the region where outside-tumor cells could have been observed: the
#' tissue outline minus the tumor polygons. Leaving tumor or empty background in
#' the window underestimates intensity and biases g(r) upward.
#'
#' @param coords_df Data frame with `x` and `y` for all cells (tumor included).
#' @param polygons sf polygon object representing tumor regions.
#' @param tissue_window Tissue outline used when `window_obj` is `NULL`:
#'   `"concave"` (concave hull of all cells), `"convex"` (convex hull), or
#'   `"bbox"` (bounding box).
#' @param window_obj Optional `spatstat.geom::owin` tissue window. Tumor polygons
#'   are subtracted from it.
#' @return A `spatstat.geom::owin` window.
#' @keywords internal
build_outside_window <- function(coords_df,
                                 polygons,
                                 tissue_window = c("concave", "convex", "bbox"),
                                 window_obj = NULL) {
  tissue_window <- match.arg(tissue_window)

  if (!is.null(window_obj)) {
    if (!inherits(window_obj, "owin")) {
      stop("`window_obj` must be a spatstat.geom::owin object.", call. = FALSE)
    }
    tissue <- sf::st_as_sfc(window_obj)
  } else {
    pts <- as.matrix(coords_df[, c("x", "y")])
    tissue <- switch(
      tissue_window,
      bbox = sf::st_as_sfc(sf::st_bbox(c(
        xmin = min(pts[, 1]), ymin = min(pts[, 2]),
        xmax = max(pts[, 1]), ymax = max(pts[, 2])
      ))),
      convex = sf::st_convex_hull(sf::st_combine(
        sf::st_as_sf(coords_df[, c("x", "y")], coords = c("x", "y"))
      )),
      concave = sf::st_sfc(sf::st_polygon(list(
        concaveman::concaveman(pts, concavity = 2)
      )))
    )
    # Hull edges pass through cells; pad slightly so those cells fall inside.
    if (tissue_window != "bbox") {
      pad <- 1e-6 * max(diff(range(pts[, 1])), diff(range(pts[, 2])), 1)
      tissue <- sf::st_buffer(tissue, pad, nQuadSegs = 1)
    }
  }

  tissue <- sf::st_set_crs(sf::st_geometry(tissue), NA)
  tumor <- sf::st_set_crs(sf::st_geometry(polygons), NA)
  tumor <- sf::st_union(sf::st_make_valid(tumor))
  outside <- sf::st_difference(sf::st_union(tissue), tumor)

  if (length(outside) == 0 || isTRUE(sf::st_is_empty(outside)) ||
      as.numeric(sf::st_area(outside)) <= 0) {
    stop("The tissue window minus the tumor polygons is empty.", call. = FALSE)
  }

  spatstat.geom::as.owin(outside)
}

# Band index of each outside cell by distance to the nearest tumor polygon.
compute_distance_bands <- function(x, y, polygons, band_width) {
  if (!is.finite(band_width)) {
    return(rep(0L, length(x)))
  }
  pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"))
  tumor <- sf::st_union(sf::st_make_valid(sf::st_set_crs(sf::st_geometry(polygons), NA)))
  d <- as.numeric(sf::st_distance(pts, tumor))
  as.integer(floor(d / band_width))
}

# Validate a pcfcross bandwidth: NULL (spatstat default) or one positive number.
check_bw <- function(bw) {
  if (!is.null(bw) && !(is.numeric(bw) && length(bw) == 1 && is.finite(bw) && bw > 0)) {
    stop("`bw` must be a single positive number or NULL.", call. = FALSE)
  }
  invisible(bw)
}

# Simulated patterns with cell-type labels permuted among cells in the same band.
simulate_band_permutations <- function(ppp_obj, bands, nsim) {
  marks0 <- spatstat.geom::marks(ppp_obj)
  idx_by_band <- split(seq_along(bands), bands)
  lapply(seq_len(nsim), function(s) {
    m <- marks0
    for (ix in idx_by_band) {
      m[ix] <- marks0[ix[sample.int(length(ix))]]
    }
    spatstat.geom::setmarks(ppp_obj, m)
  })
}

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
#' @param window_obj Optional `spatstat.geom::owin` tissue window. Tumor polygons
#'   are subtracted from it. If `NULL`, the tissue outline is built from all cell
#'   coordinates according to `tissue_window`.
#' @param tissue_window Tissue outline used when `window_obj` is `NULL`; see
#'   [build_outside_window()].
#' @return A list with `ppp_obj`, `cells_outside`, and `window`.
#' @keywords internal
build_outside_ppp <- function(seurat_obj,
                              polygons,
                              celltype_col,
                              fov_name = "fov",
                              min_cells = 10,
                              window_obj = NULL,
                              tissue_window = "concave") {
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

  window_obj <- build_outside_window(coords_df, polygons, tissue_window, window_obj)
  in_window <- spatstat.geom::inside.owin(cells_outside$x, cells_outside$y, window_obj)
  if (any(!in_window)) {
    warning(sum(!in_window), " outside-polygon cell(s) fall outside the PCF window ",
            "and were dropped.", call. = FALSE)
    cells_outside <- cells_outside[in_window, , drop = FALSE]
  }

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
#' @param distance Target distance (microns) at which to extract the PCF value.
#'   g is evaluated at exactly this distance.
#' @param min_cells Minimum number of cells required per cell type.
#' @param celltype_col Metadata column containing cell type labels.
#' @param fov_name Field of view/image name passed to `Seurat::GetTissueCoordinates`.
#' @param include_self_pairs Logical; include same-type pairs.
#' @param ordered_pairs Logical; if `TRUE`, computes ordered pairs (A,B and B,A).
#'   If `FALSE`, computes each pair once.
#' @param window_obj Optional `spatstat.geom::owin` tissue window; tumor polygons
#'   are subtracted from it.
#' @param show_progress Logical; show a text progress bar.
#' @param tissue_window Tissue outline used when `window_obj` is `NULL`:
#'   `"concave"` (default; concave hull of all cells), `"convex"`, or `"bbox"`.
#'   Tumor polygons are always removed from the window.
#' @param bw Kernel bandwidth (microns) passed to `pcfcross`. Use the same value
#'   as in [plot_pcf_outside_pair_envelopes()] so heatmap values match the
#'   envelope curves. `NULL` uses the `spatstat` default.
#' @return Data frame with columns `type1`, `type2`, `g_value`, `n_type1`, `n_type2`.
#' @export
calculate_pcf_outside_summary <- function(seurat_obj,
                                          polygons,
                                          distance = 30,
                                          min_cells = 10,
                                          celltype_col,
                                          fov_name = "fov",
                                          include_self_pairs = TRUE,
                                          ordered_pairs = FALSE,
                                          window_obj = NULL,
                                          show_progress = interactive(),
                                          tissue_window = c("concave", "convex", "bbox"),
                                          bw = 10) {
  if (!requireNamespace("spatstat.explore", quietly = TRUE)) {
    stop("Package 'spatstat.explore' is required. Install with: install.packages('spatstat.explore')")
  }
  tissue_window <- match.arg(tissue_window)
  check_bw(bw)
  if (!(is.numeric(distance) && length(distance) == 1 &&
        is.finite(distance) && distance > 0)) {
    stop("`distance` must be a single positive number.", call. = FALSE)
  }
  # pcfcross needs an evenly spaced grid from 0; ending it at `distance` gives g there exactly.
  r_grid <- seq(0, distance, length.out = 257)

  prep <- build_outside_ppp(
    seurat_obj = seurat_obj,
    polygons = polygons,
    celltype_col = celltype_col,
    fov_name = fov_name,
    min_cells = min_cells,
    window_obj = window_obj,
    tissue_window = tissue_window
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
      pcf_ij <- spatstat.explore::pcfcross(ppp_obj, i = type1, j = type2,
                                           r = r_grid, bw = bw)
      g_val <- pcf_ij$iso[length(r_grid)]

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
#' Two null models are available:
#' * `"csr"`: complete spatial randomness. Points are placed uniformly at random
#'   in the outside-tumor window with independent labels. This tests clustering of
#'   the cells themselves as well as co-localization, and tissue is rarely uniform.
#' * `"band_permutation"`: cell positions are kept and cell-type labels are
#'   permuted among cells in the same distance band from the tumor boundary
#'   (bands of width `band_width`). This tests whether two types co-localize
#'   beyond what their shared distance-to-tumor gradients explain. Use
#'   `band_width = Inf` for plain random labelling across all outside cells.
#'
#' Two envelope types are available:
#' * `"global"` (default): a simultaneous envelope (`spatstat` `global = TRUE`,
#'   maximum absolute deviation). The observed curve leaves it anywhere in
#'   `r_range` with probability `1 / (nsim + 1)` under the null, so it can be read
#'   across all distances at once without multiple testing. The band has the
#'   same width at every distance. A matching MAD test p-value, computed from the
#'   same simulations, is stored in the `p_value` attribute. With
#'   `null_model = "band_permutation"`, `2 * nsim` permutations are run: the
#'   second set estimates the null mean curve.
#' * `"pointwise"`: separate envelopes at each distance, with a two-sided
#'   significance level of `2 / (nsim + 1)` at each single r. Checking many
#'   distances inflates the false-positive rate.
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param polygons sf polygon object representing ROIs.
#' @param celltype_col Metadata column containing cell type labels.
#' @param fov_name Field of view/image name passed to `Seurat::GetTissueCoordinates`.
#' @param min_cells Minimum number of cells required per cell type.
#' @param include_self_pairs Logical; include same-type pairs.
#' @param ordered_pairs Logical; if `TRUE`, computes ordered pairs.
#' @param nsim Number of simulations for envelopes (required). Global envelopes
#'   have significance level `1 / (nsim + 1)` (0.05 for `nsim = 19`, 0.01 for
#'   `nsim = 99`); pointwise envelopes have `2 / (nsim + 1)` at each r.
#' @param bw Kernel bandwidth (microns) passed to `pcfcross`. Use the same value
#'   as in [calculate_pcf_outside_summary()] so the envelope curves match the
#'   heatmap values. `NULL` uses the `spatstat` default.
#' @param window_obj Optional `spatstat.geom::owin` tissue window; tumor polygons
#'   are subtracted from it.
#' @param file Optional base PNG file path (e.g. `"pcf_envelopes.png"`). One PNG is
#'   written per pair, with the pair name appended to the file stem.
#' @param width,height PNG size in inches.
#' @param res PNG resolution in pixels per inch.
#' @param verbose Logical; print progress messages.
#' @param tissue_window Tissue outline used when `window_obj` is `NULL`:
#'   `"concave"` (default; concave hull of all cells), `"convex"`, or `"bbox"`.
#'   Tumor polygons are always removed from the window.
#' @param null_model Null model for the envelopes: `"csr"` or
#'   `"band_permutation"`. See Details.
#' @param band_width Band width (microns) for `null_model = "band_permutation"`.
#' @param envelope_type `"global"` (default) or `"pointwise"`. See Details.
#' @param r_range Distance range `c(lower, upper)` in microns. `upper` sets the
#'   largest distance computed (default: the `spatstat` default maximum). For
#'   global envelopes the test covers `[lower, upper]`; `lower` defaults to `bw`
#'   (or 5% of `upper` when `bw` is `NULL`) because g(r) is unstable near 0.
#' @return Named list of envelope objects. Each has attributes `null_model`,
#'   `envelope_type`, `alpha` (significance level), and `r_range`; global
#'   envelopes also have `p_value` (MAD test).
#' @export
plot_pcf_outside_pair_envelopes <- function(seurat_obj,
                                            polygons,
                                            celltype_col,
                                            fov_name = "fov",
                                            min_cells = 10,
                                            include_self_pairs = FALSE,
                                            ordered_pairs = FALSE,
                                            nsim,
                                            bw = 10,
                                            window_obj = NULL,
                                            file = NULL,
                                            width = 7,
                                            height = 5,
                                            res = 150,
                                            verbose = TRUE,
                                            tissue_window = c("concave", "convex", "bbox"),
                                            null_model = c("csr", "band_permutation"),
                                            band_width = 25,
                                            envelope_type = c("global", "pointwise"),
                                            r_range = NULL) {
  if (!requireNamespace("spatstat.explore", quietly = TRUE)) {
    stop("Package 'spatstat.explore' is required. Install with: install.packages('spatstat.explore')")
  }
  tissue_window <- match.arg(tissue_window)
  null_model <- match.arg(null_model)
  envelope_type <- match.arg(envelope_type)
  if (missing(nsim)) {
    stop("`nsim` (number of envelope simulations) must be supplied, e.g. nsim = 19 ",
         "for a global envelope at significance level 0.05 (nsim = 39 for pointwise).",
         call. = FALSE)
  }
  if (!(is.numeric(nsim) && length(nsim) == 1 && is.finite(nsim) &&
        nsim >= 1 && nsim == round(nsim))) {
    stop("`nsim` must be a single positive whole number.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  check_bw(bw)
  if (null_model == "band_permutation" &&
      !(is.numeric(band_width) && length(band_width) == 1 &&
        !is.na(band_width) && band_width > 0)) {
    stop("`band_width` must be a single positive number (or Inf).", call. = FALSE)
  }
  if (!is.null(r_range) &&
      !(is.numeric(r_range) && length(r_range) == 2 && all(is.finite(r_range)) &&
        r_range[1] >= 0 && r_range[2] > r_range[1])) {
    stop("`r_range` must be NULL or c(lower, upper) with 0 <= lower < upper.", call. = FALSE)
  }
  alpha <- if (envelope_type == "global") 1 / (nsim + 1) else 2 / (nsim + 1)

  prep <- build_outside_ppp(
    seurat_obj = seurat_obj,
    polygons = polygons,
    celltype_col = celltype_col,
    fov_name = fov_name,
    min_cells = min_cells,
    window_obj = window_obj,
    tissue_window = tissue_window
  )

  ppp_obj <- prep$ppp_obj
  cell_types <- levels(spatstat.geom::marks(ppp_obj))

  # One r grid for every pair, so all curves and global tests cover the same distances.
  r_max <- if (is.null(r_range)) {
    max(spatstat.explore::pcfcross(ppp_obj, i = cell_types[1], j = cell_types[2], bw = bw)$r)
  } else {
    r_range[2]
  }
  r_grid <- seq(0, r_max, length.out = 513)
  r_lower <- if (!is.null(r_range)) r_range[1] else if (!is.null(bw)) min(bw, r_max / 2) else r_max / 20
  test_range <- c(r_lower, r_max)

  # One shared set of simulations so every pair is tested against the same null.
  # Global envelopes without a known theoretical curve need a second nsim set for the null mean.
  sims <- NULL
  if (null_model == "band_permutation") {
    bands <- compute_distance_bands(ppp_obj$x, ppp_obj$y, polygons, band_width)
    n_patterns <- if (envelope_type == "global") 2L * nsim else nsim
    sims <- simulate_band_permutations(ppp_obj, bands, n_patterns)
  }

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
    env_list[k] <- list(tryCatch({
      env_obj <- spatstat.explore::envelope(
        ppp_obj,
        fun = spatstat.explore::pcfcross,
        i = type1,
        j = type2,
        nsim = nsim,
        simulate = sims,
        r = r_grid,
        bw = bw,
        global = envelope_type == "global",
        ginterval = if (envelope_type == "global") test_range else NULL,
        savefuns = envelope_type == "global",
        verbose = FALSE
      )
      attr(env_obj, "null_model") <- null_model
      attr(env_obj, "envelope_type") <- envelope_type
      attr(env_obj, "alpha") <- alpha
      attr(env_obj, "r_range") <- test_range
      if (envelope_type == "global") {
        attr(env_obj, "p_value") <- spatstat.explore::mad.test(
          env_obj, rinterval = test_range, verbose = FALSE
        )$p.value
      }

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
    }, finally = if (png_open) grDevices::dev.off()))
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
#' env_list <- plot_pcf_outside_pair_envelopes(seurat_obj, auto_polygons,
#'                                             celltype_col = "cellType", nsim = 19)
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
  null_attr <- attr(env_obj, "null_model")
  null_label <- switch(
    if (is.null(null_attr)) "" else null_attr,
    csr = "CSR null",
    band_permutation = "band-permutation null",
    NULL
  )
  subtitle <- if (is.null(nsim)) {
    "Outside-polygon PCF"
  } else {
    sprintf("Outside-polygon PCF  |  envelope from %d simulations", as.integer(nsim))
  }
  if (!is.null(null_label)) {
    subtitle <- paste0(subtitle, " (", null_label, ")")
  }
  env_type <- attr(env_obj, "envelope_type")
  alpha <- attr(env_obj, "alpha")
  p_value <- attr(env_obj, "p_value")
  if (!is.null(env_type) && !is.null(alpha)) {
    subtitle <- paste0(subtitle, "\n", env_type, " envelope, alpha = ", signif(alpha, 2))
    if (!is.null(p_value)) {
      subtitle <- paste0(subtitle, "  |  MAD test p = ", signif(p_value, 2))
    }
  }
  env_df <- as.data.frame(env_obj)
  # A global envelope is only valid over the distance range it was computed for.
  band_df <- env_df
  test_range <- attr(env_obj, "r_range")
  if (identical(env_type, "global") && !is.null(test_range)) {
    band_df <- env_df[env_df$r >= test_range[1] & env_df$r <= test_range[2], , drop = FALSE]
  }

  ggplot2::ggplot(env_df, ggplot2::aes(x = .data$r)) +
    ggplot2::geom_ribbon(
      data = band_df,
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
