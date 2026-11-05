#' Calculate Pair Correlation Within Distance Bands Around Tumor Polygons
#'
#' Splits the outside-tumor region into bands that follow the tumor contour
#' (band 1 = `band_breaks[1]` to `band_breaks[2]` from the tumor boundary, and so
#' on) and computes the cross pair correlation function g(r) for every
#' cell-type pair separately in each band. This shows whether co-localization
#' of two cell types changes with distance from the tumor.
#'
#' Band membership uses each cell's exact distance to the nearest tumor polygon
#' (the same definition as the `"band_permutation"` null in
#' [plot_pcf_outside_pair_envelopes()]). Band shapes follow the irregular tumor
#' outline and are clipped to the tissue window.
#'
#' Two approaches are available:
#' * `"restricted"`: each band is its own observation window and only cells in
#'   that band are used, for both cell types. g(r) answers "within this band, are
#'   type-2 cells concentrated around type-1 cells, relative to the band's own
#'   cell mix?". Bands are thin, so circles of radius r larger than the band width
#'   mostly fall outside the band and estimates become noisy; keep `distance`
#'   at or below roughly the band width.
#' * `"anchored"`: type-1 ("anchor") cells come from the band; type-2
#'   ("neighbor") cells can be anywhere outside the tumor, and the window is the
#'   whole outside-tumor region. g(r) answers "around type-1 cells located in
#'   this band, how dense are type-2 cells at distance r, relative to type-2
#'   density over the whole outside region?". Not symmetric in the two types, so
#'   ordered pairs are used by default. Larger `distance` values are usable, but
#'   with `null_model = "csr"` the tumor gradient is not accounted for (see below).
#'
#' Null models (Monte Carlo, `nsim` simulations):
#' * `"label_shuffle"`: cell positions are kept and cell-type labels are randomly
#'   permuted. With `"restricted"`, labels are permuted among cells in the band.
#'   With `"anchored"`, labels are permuted within every band (cells beyond the
#'   last band form their own group), so each band keeps its cell-type
#'   composition and the distance-to-tumor gradient is preserved. Tests
#'   co-localization beyond what the band's composition and cell layout explain.
#' * `"csr"`: complete spatial randomness. With `"restricted"`, every cell in the
#'   band is moved to a uniformly random position within the band, keeping its
#'   type. With `"anchored"`, anchor cells keep their real positions and every
#'   other outside-tumor cell is moved to a uniformly random position in the
#'   outside window. CSR also flags ordinary clumping of cells, and in
#'   `"anchored"` mode it also flags any distance-to-tumor gradient.
#' * `"none"`: observed g only.
#'
#' The p-value is a two-sided Monte Carlo test at each `distance`:
#' `(1 + k) / (nsim_ok + 1)`, where `k` is the number of simulations with
#' `abs(g_sim - null_mean) >= abs(g_obs - null_mean)` and `nsim_ok` is the number
#' of simulations giving a finite g. `p_adj`
#' adjusts these across all rows of the result (all bands, pairs and distances)
#' with `p_adjust_method`.
#'
#' @param seurat_obj Seurat object with spatial coordinates.
#' @param polygons sf polygon object representing tumor regions.
#' @param celltype_col Metadata column containing cell type labels.
#' @param band_width Band width in coordinate units (microns). Ignored when
#'   `band_breaks` is supplied.
#' @param max_distance Outer edge of the last band. `NULL` covers all outside
#'   cells (rounded up to a whole band). Ignored when `band_breaks` is supplied.
#' @param band_breaks Optional increasing numeric vector of band edges, e.g.
#'   `c(0, 25, 50, 100, 200)`. Overrides `band_width` and `max_distance`. Cells
#'   closer than the first or at least the last break are in no band.
#' @param approach `"restricted"` or `"anchored"`. See Details.
#' @param distance Cell-to-cell distance(s) r at which g is reported. One value
#'   or a vector.
#' @param null_model `"label_shuffle"`, `"csr"`, or `"none"`. See Details.
#' @param nsim Number of simulations for the null model.
#' @param bw Kernel bandwidth passed to `spatstat.explore::pcfmulti()`. `NULL`
#'   uses the `spatstat` default.
#' @param min_cells Minimum cells of a type needed to compute g for it: per band
#'   for both types with `"restricted"`; per band for the anchor type and over
#'   the whole outside region for the neighbor type with `"anchored"`. Types
#'   with fewer than `min_cells` outside cells overall are dropped entirely.
#' @param fov_name Field of view/image name passed to `Seurat::GetTissueCoordinates`.
#' @param include_self_pairs Logical; include same-type pairs.
#' @param ordered_pairs Logical; compute both A-B and B-A. `NULL` (default) uses
#'   `TRUE` for `"anchored"` and `FALSE` for `"restricted"`.
#' @param window_obj Optional `spatstat.geom::owin` tissue window; tumor polygons
#'   are subtracted from it.
#' @param tissue_window Tissue outline used when `window_obj` is `NULL`:
#'   `"concave"`, `"convex"`, or `"bbox"`. See [build_outside_window()].
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param show_progress Logical; show a text progress bar.
#' @return Data frame with one row per band, pair and distance: `band`,
#'   `band_inner`, `band_outer`, `band_label`, `type1`, `type2`, `distance`,
#'   `g_value`, `n_type1`, `n_type2` (cells used: in the band, or for the
#'   anchored neighbor type, in the whole outside region), and when a null model
#'   is used `null_mean`, `null_lo`, `null_hi` (2.5% and 97.5% simulation
#'   quantiles), `p_value`, `p_adj`. Attributes: `approach`, `null_model`,
#'   `nsim`, `bw`, `band_breaks`, `ordered_pairs`, and `bands` (an sf object of
#'   band polygons for plotting).
#' @export
#'
#' @examples
#' \dontrun{
#' res <- calculate_pcf_by_band(seurat_obj, auto_polygons, celltype_col = "cellType",
#'                              band_width = 50, max_distance = 200,
#'                              approach = "restricted", distance = 20,
#'                              null_model = "label_shuffle", nsim = 99)
#' plot_pcf_band_heatmap(res)
#' plot_pcf_band_profile(res)
#' }
calculate_pcf_by_band <- function(seurat_obj,
                                  polygons,
                                  celltype_col,
                                  band_width = 25,
                                  max_distance = NULL,
                                  band_breaks = NULL,
                                  approach = c("restricted", "anchored"),
                                  distance = 20,
                                  null_model = c("label_shuffle", "csr", "none"),
                                  nsim = 99,
                                  bw = 10,
                                  min_cells = 10,
                                  fov_name = "fov",
                                  include_self_pairs = TRUE,
                                  ordered_pairs = NULL,
                                  window_obj = NULL,
                                  tissue_window = c("concave", "convex", "bbox"),
                                  p_adjust_method = "BH",
                                  show_progress = interactive()) {
  for (p in c("spatstat.explore", "spatstat.geom")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required. Install with: install.packages('", p, "')",
           call. = FALSE)
    }
  }
  approach <- match.arg(approach)
  null_model <- match.arg(null_model)
  tissue_window <- match.arg(tissue_window)
  if (null_model == "csr" && !requireNamespace("spatstat.random", quietly = TRUE)) {
    stop("Package 'spatstat.random' is required for null_model = 'csr'. ",
         "Install with: install.packages('spatstat.random')", call. = FALSE)
  }
  check_bw(bw)
  if (!(is.numeric(distance) && length(distance) >= 1 && all(is.finite(distance)) &&
        all(distance > 0))) {
    stop("`distance` must be one or more positive numbers.", call. = FALSE)
  }
  distance <- sort(unique(distance))
  if (null_model != "none" &&
      !(is.numeric(nsim) && length(nsim) == 1 && is.finite(nsim) &&
        nsim >= 1 && nsim == round(nsim))) {
    stop("`nsim` must be a single positive whole number.", call. = FALSE)
  }
  nsim <- if (null_model == "none") 0L else as.integer(nsim)
  if (is.null(ordered_pairs)) {
    ordered_pairs <- approach == "anchored"
  }

  prep <- build_outside_ppp(
    seurat_obj = seurat_obj,
    polygons = polygons,
    celltype_col = celltype_col,
    fov_name = fov_name,
    min_cells = min_cells,
    window_obj = window_obj,
    tissue_window = tissue_window
  )
  X <- prep$ppp_obj
  W <- prep$window
  marks_all <- spatstat.geom::marks(X)
  cell_types <- levels(marks_all)

  tumor <- sf::st_union(sf::st_make_valid(sf::st_set_crs(sf::st_geometry(polygons), NA)))
  pts <- sf::st_as_sf(data.frame(x = X$x, y = X$y), coords = c("x", "y"))
  d_tumor <- as.numeric(sf::st_distance(pts, tumor))

  band_breaks <- resolve_band_breaks(band_breaks, band_width, max_distance, d_tumor)
  n_bands <- length(band_breaks) - 1L
  band_idx <- findInterval(d_tumor, band_breaks, left.open = FALSE, rightmost.closed = FALSE)
  band_idx[band_idx < 1L | band_idx > n_bands] <- NA_integer_

  if (approach == "restricted" && max(distance) > min(diff(band_breaks))) {
    warning("With approach = 'restricted', g at distances larger than the band width ",
            "(", signif(min(diff(band_breaks)), 3), ") relies heavily on edge ",
            "correction and is unreliable.", call. = FALSE)
  }

  outside_sf <- sf::st_set_crs(sf::st_as_sfc(W), NA)
  band_geoms <- build_band_geometries(tumor, outside_sf, band_breaks)
  band_labels <- paste0(format_num(band_breaks[-length(band_breaks)]), "\u2013",
                        format_num(band_breaks[-1]))

  pair_grid <- make_pcf_pair_grid(cell_types, include_self_pairs, ordered_pairs)
  r_grid <- seq(0, max(distance), length.out = 257)

  g_at <- function(Xp, I, J) {
    if (sum(I) < 1 || sum(J) < 1 || (sum(I) == 1 && sum(J) == 1 && any(I & J))) {
      return(rep(NA_real_, length(distance)))
    }
    out <- tryCatch({
      fv <- spatstat.explore::pcfmulti(Xp, I, J, r = r_grid, bw = bw,
                                       correction = "Ripley")
      stats::approx(fv$r, fv$iso, xout = distance, rule = 2)$y
    }, error = function(e) rep(NA_real_, length(distance)))
    out[!is.finite(out)] <- NA_real_
    out
  }

  # Stratum per outside cell for anchored label shuffling (cells beyond the last band = 0).
  strata_all <- ifelse(is.na(band_idx), 0L, band_idx)

  rows <- list()
  n_dropped <- 0L
  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 0, max = n_bands, style = 3)
    on.exit({
      close(pb)
      cat("\n")
    }, add = TRUE)
  }

  # Label-shuffled marks for the anchored approach are shared across bands.
  if (approach == "anchored" && null_model == "label_shuffle") {
    shuffled_marks <- lapply(seq_len(nsim), function(s) permute_within(marks_all, strata_all))
  }

  for (k in seq_len(n_bands)) {
    if (!is.null(pb)) utils::setTxtProgressBar(pb, k)

    Wk <- tryCatch(
      if (is.null(band_geoms[[k]])) NULL else spatstat.geom::as.owin(band_geoms[[k]]),
      error = function(e) NULL
    )

    if (approach == "restricted") {
      in_k <- which(!is.na(band_idx) & band_idx == k)
      if (!is.null(Wk) && length(in_k) > 0) {
        # Buffered band edges approximate arcs with straight segments.
        ok <- spatstat.geom::inside.owin(X$x[in_k], X$y[in_k], Wk)
        n_dropped <- n_dropped + sum(!ok)
        in_k <- in_k[ok]
      }
      Xk <- if (is.null(Wk) || length(in_k) == 0) NULL else
        spatstat.geom::ppp(X$x[in_k], X$y[in_k], window = Wk,
                           marks = marks_all[in_k], checkdup = FALSE)
      mk <- if (is.null(Xk)) factor(character(0), levels = cell_types) else
        spatstat.geom::marks(Xk)
      n_k <- table(factor(mk, levels = cell_types))

      sims_k <- NULL
      if (!is.null(Xk) && nsim > 0) {
        sims_k <- lapply(seq_len(nsim), function(s) {
          if (null_model == "label_shuffle") {
            spatstat.geom::setmarks(Xk, mk[sample.int(length(mk))])
          } else {
            Y <- spatstat.random::runifpoint(Xk$n, win = Wk)
            spatstat.geom::setmarks(Y, mk)
          }
        })
      }

      for (p in seq_len(nrow(pair_grid))) {
        t1 <- pair_grid$type1[p]
        t2 <- pair_grid$type2[p]
        enough <- !is.null(Xk) && n_k[[t1]] >= min_cells && n_k[[t2]] >= min_cells
        g_obs <- if (enough) g_at(Xk, mk == t1, mk == t2) else rep(NA_real_, length(distance))
        g_sim <- if (enough && nsim > 0) {
          do.call(rbind, lapply(sims_k, function(Y) {
            my <- spatstat.geom::marks(Y)
            g_at(Y, my == t1, my == t2)
          }))
        } else NULL
        rows[[length(rows) + 1L]] <- band_rows(k, band_breaks, band_labels, t1, t2, distance,
                                               g_obs, n_k[[t1]], n_k[[t2]], g_sim, nsim)
      }
    } else {
      n_total <- table(factor(marks_all, levels = cell_types))
      in_band <- !is.na(band_idx) & band_idx == k
      for (t1 in unique(pair_grid$type1)) {
        I <- in_band & marks_all == t1
        n_anchor <- sum(I)
        pairs_t1 <- pair_grid[pair_grid$type1 == t1, , drop = FALSE]
        enough_anchor <- n_anchor >= min_cells

        csr_sims <- NULL
        if (enough_anchor && nsim > 0 && null_model == "csr") {
          csr_sims <- lapply(seq_len(nsim), function(s) {
            Y <- spatstat.random::runifpoint(sum(!I), win = W)
            xy <- data.frame(x = X$x, y = X$y)
            xy$x[!I] <- Y$x
            xy$y[!I] <- Y$y
            spatstat.geom::ppp(xy$x, xy$y, window = W, marks = marks_all, checkdup = FALSE)
          })
        }

        for (p in seq_len(nrow(pairs_t1))) {
          t2 <- pairs_t1$type2[p]
          enough <- enough_anchor && n_total[[t2]] >= min_cells
          g_obs <- if (enough) g_at(X, I, marks_all == t2) else rep(NA_real_, length(distance))
          g_sim <- NULL
          if (enough && nsim > 0) {
            g_sim <- if (null_model == "label_shuffle") {
              do.call(rbind, lapply(shuffled_marks, function(ms) {
                g_at(X, in_band & ms == t1, ms == t2)
              }))
            } else {
              do.call(rbind, lapply(csr_sims, function(Y) g_at(Y, I, marks_all == t2)))
            }
          }
          rows[[length(rows) + 1L]] <- band_rows(k, band_breaks, band_labels, t1, t2, distance,
                                                 g_obs, n_anchor, n_total[[t2]], g_sim, nsim)
        }
      }
    }
  }

  if (n_dropped > 0) {
    warning(n_dropped, " cell(s) on a band edge fell just outside the band polygon ",
            "and were dropped.", call. = FALSE)
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- band_rows(integer(0), band_breaks, band_labels, character(0), character(0),
                     numeric(0), numeric(0), integer(0), integer(0), NULL, nsim)
  }
  if (nsim > 0) {
    out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust_method)
  }
  out$band_label <- factor(out$band_label, levels = band_labels)
  rownames(out) <- NULL

  bands_sf <- sf::st_sf(
    band = seq_len(n_bands),
    band_label = factor(band_labels, levels = band_labels),
    geometry = sf::st_sfc(lapply(band_geoms, function(g) {
      if (is.null(g)) sf::st_polygon() else sf::st_geometry(g)[[1]]
    }))
  )

  attr(out, "approach") <- approach
  attr(out, "null_model") <- null_model
  attr(out, "nsim") <- nsim
  attr(out, "bw") <- bw
  attr(out, "band_breaks") <- band_breaks
  attr(out, "ordered_pairs") <- ordered_pairs
  attr(out, "bands") <- bands_sf
  out
}

#' Plot Per-Band PCF Heatmaps
#'
#' Draws one cell-type by cell-type heatmap of g(r) per distance band from the
#' output of [calculate_pcf_by_band()].
#'
#' @param band_summary Output of [calculate_pcf_by_band()].
#' @param distance Which `distance` to plot when several were computed. `NULL`
#'   uses the smallest.
#' @param value `"g_value"` (observed g) or `"g_over_null"` (observed g divided
#'   by the null mean). `"g_over_null"` removes effects already present in the
#'   null, e.g. the tumor gradient under `"label_shuffle"` with `"anchored"`.
#'   Requires a null model.
#' @param midpoint Midpoint of the diverging color scale.
#' @param alpha Significance level for marking cells with `*` (uses `p_adj`).
#'   Set `NULL` to hide markers.
#' @return ggplot object, faceted by band.
#' @export
plot_pcf_band_heatmap <- function(band_summary, distance = NULL,
                                  value = c("g_value", "g_over_null"),
                                  midpoint = 1, alpha = 0.05) {
  value <- match.arg(value)
  df <- select_band_distance(band_summary, distance)
  if (value == "g_over_null") {
    if (!"null_mean" %in% names(df)) {
      stop("value = 'g_over_null' needs results computed with a null model.", call. = FALSE)
    }
    df$g_value <- df$g_value / df$null_mean
  }
  ordered <- isTRUE(attr(band_summary, "ordered_pairs"))
  anchored <- identical(attr(band_summary, "approach"), "anchored")

  if (!ordered) {
    flip <- df[df$type1 != df$type2, , drop = FALSE]
    flip[, c("type1", "type2")] <- flip[, c("type2", "type1")]
    flip[, c("n_type1", "n_type2")] <- flip[, c("n_type2", "n_type1")]
    df <- rbind(df, flip)
  }
  df$sig <- ""
  if (!is.null(alpha) && "p_adj" %in% names(df)) {
    df$sig <- ifelse(!is.na(df$p_adj) & df$p_adj < alpha, "*", "")
  }
  types <- sort(unique(c(df$type1, df$type2)))
  df$type1 <- factor(df$type1, levels = types)
  df$type2 <- factor(df$type2, levels = types)

  .data <- NULL
  ggplot2::ggplot(df, ggplot2::aes(x = .data$type2, y = .data$type1, fill = .data$g_value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig), size = 5, vjust = 0.75) +
    ggplot2::facet_wrap(~band_label, labeller = ggplot2::labeller(
      band_label = function(x) paste("Band", x)
    )) +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                  midpoint = midpoint, na.value = "grey90") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      title = "Pair correlation by distance band from tumor",
      subtitle = band_plot_subtitle(band_summary, unique(df$distance), alpha),
      x = if (anchored) "Neighbor cell type (anywhere outside tumor)" else "Cell type",
      y = if (anchored) "Anchor cell type (in band)" else "Cell type",
      fill = if (value == "g_over_null") "g / null mean" else "g(r)"
    )
}

#' Plot PCF Across Distance Bands
#'
#' Plots g(r) against distance from the tumor (band midpoint) for each cell-type
#' pair, with the null 95% simulation interval when available.
#'
#' @param band_summary Output of [calculate_pcf_by_band()].
#' @param pairs Optional character vector of pairs to show, written
#'   `"TypeA vs TypeB"`. `NULL` shows all pairs.
#' @param distance Which `distance` to plot when several were computed. `NULL`
#'   uses the smallest.
#' @param alpha Significance level for filled points (uses `p_adj`).
#' @return ggplot object, faceted by pair.
#' @export
plot_pcf_band_profile <- function(band_summary, pairs = NULL, distance = NULL, alpha = 0.05) {
  df <- select_band_distance(band_summary, distance)
  df$pair <- paste(df$type1, "vs", df$type2)
  if (!is.null(pairs)) {
    missing_pairs <- setdiff(pairs, df$pair)
    if (length(missing_pairs) > 0) {
      stop("Pair(s) not found: ", paste(missing_pairs, collapse = ", "),
           ". Available: ", paste(unique(df$pair), collapse = ", "), call. = FALSE)
    }
    df <- df[df$pair %in% pairs, , drop = FALSE]
  }
  df$band_mid <- (df$band_inner + df$band_outer) / 2
  df$significant <- if ("p_adj" %in% names(df)) !is.na(df$p_adj) & df$p_adj < alpha else FALSE

  .data <- NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$band_mid, y = .data$g_value))
  if (all(c("null_lo", "null_hi") %in% names(df))) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$null_lo, ymax = .data$null_hi),
                                  fill = "grey70", alpha = 0.5)
  }
  p +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$significant), color = "steelblue",
                        size = 2.5, fill = "steelblue") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19),
                                labels = c(`FALSE` = "not significant",
                                           `TRUE` = paste0("p_adj < ", alpha)),
                                drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~pair, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Pair correlation vs. distance from tumor",
      subtitle = band_plot_subtitle(band_summary, unique(df$distance), alpha),
      x = "Distance from tumor boundary (band midpoint)",
      y = "g(r)"
    )
}

# ---- internal helpers -------------------------------------------------------

resolve_band_breaks <- function(band_breaks, band_width, max_distance, d_tumor) {
  if (!is.null(band_breaks)) {
    if (!(is.numeric(band_breaks) && length(band_breaks) >= 2 &&
          all(is.finite(band_breaks)) && all(band_breaks >= 0) &&
          all(diff(band_breaks) > 0))) {
      stop("`band_breaks` must be an increasing vector of at least two non-negative numbers.",
           call. = FALSE)
    }
    return(band_breaks)
  }
  if (!(is.numeric(band_width) && length(band_width) == 1 && is.finite(band_width) &&
        band_width > 0)) {
    stop("`band_width` must be a single positive number.", call. = FALSE)
  }
  if (is.null(max_distance)) {
    # Bands are half-open, so add one if the farthest cell sits exactly on the last edge.
    n <- max(1L, ceiling(max(d_tumor) / band_width))
    if (n * band_width <= max(d_tumor)) n <- n + 1L
  } else if (!(is.numeric(max_distance) && length(max_distance) == 1 &&
               is.finite(max_distance) && max_distance > 0)) {
    stop("`max_distance` must be NULL or a single positive number.", call. = FALSE)
  } else {
    n <- max(1L, ceiling(max_distance / band_width - 1e-9))
  }
  seq(0, n * band_width, by = band_width)
}

build_band_geometries <- function(tumor, outside_sf, band_breaks) {
  lapply(seq_len(length(band_breaks) - 1L), function(k) {
    outer <- sf::st_buffer(tumor, band_breaks[k + 1], nQuadSegs = 32)
    ring <- if (band_breaks[k] > 0) {
      sf::st_difference(outer, sf::st_buffer(tumor, band_breaks[k], nQuadSegs = 32))
    } else {
      outer
    }
    g <- sf::st_intersection(sf::st_union(outside_sf), ring)
    if (length(g) == 0 || all(sf::st_is_empty(g))) return(NULL)
    g <- sf::st_make_valid(g)
    if (any(sf::st_geometry_type(g) == "GEOMETRYCOLLECTION")) {
      g <- sf::st_collection_extract(g, "POLYGON")
    }
    g <- sf::st_union(g)
    if (length(g) == 0 || all(sf::st_is_empty(g)) || as.numeric(sf::st_area(g)) <= 0) NULL else g
  })
}

make_pcf_pair_grid <- function(cell_types, include_self_pairs, ordered_pairs) {
  if (ordered_pairs) {
    grid <- expand.grid(type1 = cell_types, type2 = cell_types, stringsAsFactors = FALSE)
    if (!include_self_pairs) grid <- grid[grid$type1 != grid$type2, , drop = FALSE]
    return(grid[order(grid$type1, grid$type2), , drop = FALSE])
  }
  grid <- if (length(cell_types) >= 2) {
    pr <- utils::combn(cell_types, 2)
    data.frame(type1 = pr[1, ], type2 = pr[2, ], stringsAsFactors = FALSE)
  } else {
    data.frame(type1 = character(0), type2 = character(0), stringsAsFactors = FALSE)
  }
  if (include_self_pairs) {
    grid <- rbind(data.frame(type1 = cell_types, type2 = cell_types,
                             stringsAsFactors = FALSE), grid)
  }
  grid
}

permute_within <- function(marks, strata) {
  out <- marks
  for (ix in split(seq_along(strata), strata)) {
    out[ix] <- marks[ix[sample.int(length(ix))]]
  }
  out
}

band_rows <- function(k, band_breaks, band_labels, t1, t2, distance, g_obs,
                      n1, n2, g_sim, nsim) {
  nd <- length(distance)
  df <- data.frame(
    band = rep(k, length.out = nd),
    band_inner = if (length(k)) band_breaks[k] else numeric(0),
    band_outer = if (length(k)) band_breaks[k + 1] else numeric(0),
    band_label = if (length(k)) band_labels[k] else character(0),
    type1 = rep(t1, length.out = nd),
    type2 = rep(t2, length.out = nd),
    distance = distance,
    g_value = g_obs,
    n_type1 = rep(as.integer(n1), length.out = nd),
    n_type2 = rep(as.integer(n2), length.out = nd),
    stringsAsFactors = FALSE
  )
  if (nsim > 0) {
    df$null_mean <- NA_real_
    df$null_lo <- NA_real_
    df$null_hi <- NA_real_
    df$p_value <- NA_real_
    if (!is.null(g_sim) && nd > 0) {
      g_sim <- matrix(g_sim, ncol = nd)
      for (m in seq_len(nd)) {
        s <- g_sim[, m]
        s <- s[is.finite(s)]
        if (length(s) == 0) next
        mu <- mean(s)
        df$null_mean[m] <- mu
        df$null_lo[m] <- unname(stats::quantile(s, 0.025))
        df$null_hi[m] <- unname(stats::quantile(s, 0.975))
        if (is.finite(g_obs[m])) {
          df$p_value[m] <- (1 + sum(abs(s - mu) >= abs(g_obs[m] - mu) - 1e-12)) /
            (length(s) + 1)
        }
      }
    }
  }
  df
}

select_band_distance <- function(band_summary, distance) {
  if (!is.data.frame(band_summary) || !all(c("band_label", "g_value", "distance") %in%
                                           names(band_summary))) {
    stop("`band_summary` must be the output of calculate_pcf_by_band().", call. = FALSE)
  }
  if (nrow(band_summary) == 0) stop("`band_summary` is empty.", call. = FALSE)
  avail <- sort(unique(band_summary$distance))
  if (is.null(distance)) {
    distance <- avail[1]
    if (length(avail) > 1) {
      message("Several distances computed; plotting r = ", distance,
              ". Use `distance` to choose another.")
    }
  } else if (!any(abs(avail - distance) < 1e-9)) {
    stop("`distance` ", distance, " not in results. Available: ",
         paste(avail, collapse = ", "), call. = FALSE)
  }
  band_summary[abs(band_summary$distance - distance) < 1e-9, , drop = FALSE]
}

band_plot_subtitle <- function(band_summary, distance, alpha) {
  null_model <- attr(band_summary, "null_model")
  s <- paste0(attr(band_summary, "approach"), " approach  |  r = ", format_num(distance))
  if (!is.null(null_model) && null_model != "none") {
    s <- paste0(s, "  |  null: ", null_model, " (", attr(band_summary, "nsim"), " sims)")
    if (!is.null(alpha)) s <- paste0(s, ", marked: p_adj < ", alpha)
  }
  s
}

format_num <- function(x) format(x, trim = TRUE, drop0trailing = TRUE)
