#' Default brca Atlas Cell Type Color Palette
#'
#' Returns a named vector of colors for singleR.predicted.id.brcaAtlas cell types.
#' Used internally when celltype_col is "singleR.predicted.id.brcaAtlas" and no
#' custom colors are provided.
#'
#' @return Named character vector of colors
#' @keywords internal
brca_atlas_colors <- function() {
  stats::setNames(
    c("#9DB0F5", "#676EA3",
      "#d95f0e", "orange",
      "mediumpurple1", "#a872cc", "#51047A", "#caa8e6", "#873bb3",
      "deepskyblue2",
      "violetred3",
      "#bae4bc",
      "springgreen",
      "cyan", "turquoise3", "turquoise", "#237874",
      "#BFAF00",
      "#2CA02C",
      "gold",
      "#046E29",
      "#ffffb2",
      "maroon1", "#fbb4b9",
      "dodgerblue", "mediumblue",
      "grey75",
      "#a50f15", "#F23333",
      "bisque2"),
    c("B cells Memory", "B cells Naive",
      "CAFs MSC iCAF-like", "CAFs myCAF-like",
      "Cancer Basal SC", "Cancer Cycling", "Cancer Her2 SC", "Cancer LumA SC", "Cancer LumB SC",
      "Cycling PVL",
      "Cycling T-cells",
      "Cycling_Myeloid",
      "DCs",
      "Endothelial ACKR1", "Endothelial CXCL12", "Endothelial Lymphatic LYVE1", "Endothelial RGS5",
      "Luminal Progenitors",
      "Macrophage",
      "Mature Luminal",
      "Monocyte",
      "Myoepithelial",
      "NK cells", "NKT cells",
      "PVL Differentiated", "PVL Immature",
      "Plasmablasts",
      "T cells CD4+", "T cells CD8+",
      "Unknown")
  )
}

#' Plot Distance Profile Comparison
#'
#' Plots transcriptomic distance for different inside modes
#'
#' @param distance_profile_combined Combined data frame from both modes
#' @return ggplot object
#' @importFrom magrittr %>%
#' @export
plot_distance_comparison <- function(distance_profile_combined) {
  ggplot2::ggplot(distance_profile_combined, 
                  ggplot2::aes(x = distance, y = transcriptomic_distance, 
                               color = inside_mode, linetype = inside_mode)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c("all" = "steelblue", "distance" = "darkorange"),
                                labels = c("all" = "All inside vs boundary outside",
                                           "distance" = "Boundary inside vs boundary outside")) +
    ggplot2::scale_linetype_manual(values = c("all" = "solid", "distance" = "dashed"),
                                   labels = c("all" = "All inside vs boundary outside",
                                              "distance" = "Boundary inside vs boundary outside")) +
    ggplot2::labs(
      title = "Transcriptomic Distance: Comparison of Inside Modes",
      subtitle = "How does the definition of 'inside' affect separation metrics?",
      x = "Distance from boundary (μm)",
      y = "Transcriptomic distance (1 - correlation)",
      color = "Inside definition",
      linetype = "Inside definition"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot All Metrics Faceted
#'
#' Creates faceted plot of all separation metrics
#'
#' @param distance_profile_combined Combined data frame from both modes
#' @return ggplot object
#' @importFrom magrittr %>%
#' @export
plot_all_metrics <- function(distance_profile_combined) {
  distance_profile_long <- distance_profile_combined %>%
    tidyr::pivot_longer(cols = c(transcriptomic_distance, composition_difference, de_strength),
                        names_to = "metric",
                        values_to = "value")
  
  ggplot2::ggplot(distance_profile_long, 
                  ggplot2::aes(x = distance, y = value, 
                               color = inside_mode, linetype = inside_mode)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(values = c("all" = "steelblue", "distance" = "darkorange"),
                                labels = c("all" = "All inside",
                                           "distance" = "Boundary inside")) +
    ggplot2::scale_linetype_manual(values = c("all" = "solid", "distance" = "dashed"),
                                   labels = c("all" = "All inside",
                                              "distance" = "Boundary inside")) +
    ggplot2::labs(
      title = "Genomic Separation: All Metrics Comparison",
      subtitle = "Effect of inside definition on separation metrics",
      x = "Distance from boundary (μm)",
      y = "Separation metric value",
      color = "Inside definition",
      linetype = "Inside definition"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot Cell Type Contributions Heatmap
#'
#' Shows which cell types drive differences at each distance
#'
#' @param celltype_contributions Data frame from calculate_distance_profile
#' @return ggplot object
#' @export
plot_celltype_heatmap <- function(celltype_contributions) {
  ggplot2::ggplot(celltype_contributions, 
                  ggplot2::aes(x = distance, y = celltype, fill = contribution)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "white", high = "red", midpoint = 0) +
    ggplot2::labs(
      title = "Cell Type Contribution to Composition Difference",
      subtitle = "Which cell types drive inside/outside differences at each distance?",
      x = "Distance from boundary (μm)",
      y = "Cell Type",
      fill = "Squared\ndifference"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

#' Plot Cell Type Proportions
#'
#' Shows proportion trajectories for top cell types
#'
#' @param celltype_contributions Data frame from calculate_distance_profile
#' @param top_n Number of top cell types to show
#' @param span Loess smoothing span
#' @param level Confidence level for the loess ribbon
#' @param ncol Number of facet columns. If `NULL` (default), chosen automatically
#'   by `ggplot2::facet_wrap()`.
#' @return ggplot object
#' @importFrom magrittr %>%
#' @export
plot_celltype_proportions <- function(celltype_contributions, top_n = 5,
                                      span = 0.3, level = 0.95, ncol = NULL) {
    total_by_distance <- celltype_contributions %>%
        dplyr::group_by(distance) %>%
        dplyr::summarise(total_contrib = sum(contribution), .groups = "drop")
    peak_distance <- total_by_distance$distance[which.max(total_by_distance$total_contrib)]
    
    top_celltypes <- celltype_contributions %>%
        dplyr::filter(distance == peak_distance) %>%
        dplyr::arrange(dplyr::desc(contribution)) %>%
        dplyr::slice(1:top_n) %>%
        dplyr::pull(celltype)
    
    celltype_props_long <- celltype_contributions %>%
        dplyr::filter(celltype %in% top_celltypes) %>%
        tidyr::pivot_longer(cols = c(inside_prop, outside_prop),
                            names_to = "region", values_to = "proportion") %>%
        dplyr::mutate(region = ifelse(region == "inside_prop", "Inside", "Outside"),
                      signed_distance = ifelse(region == "Inside", -distance, distance)) %>%
        dplyr::arrange(celltype, signed_distance)
    
    max_dist <- max(abs(celltype_props_long$signed_distance))
    
    ggplot2::ggplot(celltype_props_long,
                    ggplot2::aes(x = signed_distance, y = proportion)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::geom_smooth(ggplot2::aes(group = 1),
                             method = "loess", formula = y ~ x, span = span,
                             se = TRUE, level = level,
                             color = "black", fill = "grey70", linewidth = 1) +
        ggplot2::geom_point(ggplot2::aes(color = region),
                            size = 1.8, alpha = 0.9, show.legend = FALSE) +
        ggplot2::facet_wrap(~celltype, scales = "free_y", ncol = ncol) +
        ggplot2::scale_color_manual(values = c(Inside = "coral", Outside = "#A7C7E7")) +
        ggplot2::scale_x_continuous(limits = c(-max_dist, max_dist)) +
        ggplot2::labs(title = "Cell Type Proportions: Inside vs Outside",
                      subtitle = paste("Top", top_n, "most different cell types"),
                      x = "Distance from boundary (\u03bcm)\n\u2190 inside tumor | outside tumor \u2192",
                      y = "Proportion") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.title.x = ggplot2::element_text(lineheight = 1.3))
}

#' QC Visualization for Xenium Spatial Data
#'
#' Generates ROI and cell type overlays from spatial coordinates and metadata.
#'
#' @param obj Seurat object with spatial coordinates
#' @param roi_col Metadata column indicating ROI membership
#' @param celltype_col Metadata column containing cell type annotation
#' @param celltype_colors Optional named vector mapping cell types to colors
#' @param strict_celltype_colors Logical; error if present cell types are missing colors
#' @param image Optional image/FOV name
#' @param coords_df Optional pre-built coordinate data.frame with `cell_id`
#' @param point_size Point size for scatter layers
#' @param alpha Point alpha
#' @param flip_y Logical; reverse y-axis
#' @param show One of "both", "roi", "celltype"
#' @param roi_colors Named vector of colors for `FALSE` and `TRUE` ROI values
#' @param drop_na_celltype Logical; drop rows with missing cell type
#' @param grey_outside_roi Logical; draw outside ROI cells in gray on celltype panel
#' @param zoom Logical; zoom to ROI bounding box or manual limits
#' @param zoom_pad Fractional padding around ROI bounding box
#' @param zoom_fixed Logical; apply same zoom limits to both panels
#' @param zoom_xlim Optional numeric vector `c(xmin, xmax)`
#' @param zoom_ylim Optional numeric vector `c(ymin, ymax)`
#' @param legend_position Legend position passed to ggplot2 theme
#' @param legend_point_size Point size in legend keys
#' @param legend_title_size Legend title text size
#' @param legend_text_size Legend text size
#' @param legend_nrow Optional legend rows
#' @param legend_ncol Optional legend columns
#' @param collect_legends Logical; collect legends when `show = "both"`
#' @return A list containing processed data and one or more ggplot objects
#' @export
qc_plot_xenium <- function(
    obj,
    roi_col = "inside_polygon",
    celltype_col = "singleR.predicted.id.brcaAtlas",
    celltype_colors = NULL,
    strict_celltype_colors = FALSE,
    image = NULL,
    coords_df = NULL,
    point_size = 0.25,
    alpha = 0.8,
    flip_y = FALSE,
    show = c("both", "roi", "celltype"),
    roi_colors = c(`FALSE` = "grey80", `TRUE` = "red"),
    drop_na_celltype = FALSE,
    grey_outside_roi = TRUE,
    zoom = FALSE,
    zoom_pad = 0.05,
    zoom_fixed = TRUE,
    zoom_xlim = NULL,
    zoom_ylim = NULL,
    legend_position = "bottom",
    legend_point_size = 2.5,
    legend_title_size = 9,
    legend_text_size = 8,
    legend_nrow = NULL,
    legend_ncol = NULL,
    collect_legends = FALSE
) {
  show <- match.arg(show)

  # Auto-apply brca atlas color palette

  if (is.null(celltype_colors) && celltype_col == "singleR.predicted.id.brcaAtlas") {
    celltype_colors <- brca_atlas_colors()
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }
  if (show == "both" && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Please install patchwork for side-by-side plots: install.packages('patchwork')")
  }

  if (!"cell_id" %in% colnames(obj@meta.data)) {
    obj$cell_id <- colnames(obj)
  }

  df <- NULL
  if (!is.null(coords_df)) {
    df <- coords_df
  } else if (!is.null(obj@misc$tissue_coords)) {
    df <- obj@misc$tissue_coords
  } else {
    imgs <- names(obj@images)
    if (length(imgs) == 0) {
      stop("No obj@images found and no coords_df/obj@misc$tissue_coords provided.")
    }
    if (is.null(image)) {
      image <- imgs[1]
    }
    if (!image %in% imgs) {
      stop("Requested image not found. Available images: ", paste(imgs, collapse = ", "))
    }

    df <- get_cell_coords(obj, image)

    obj@misc$xenium_image_name <- image
  }

  if (!"cell_id" %in% colnames(df)) {
    stop("Coordinate df must contain a 'cell_id' column (or be buildable into one).")
  }
  if (!all(c("x", "y") %in% colnames(df))) {
    if (all(c("col", "row") %in% colnames(df))) {
      df$x <- df$col
      df$y <- df$row
    } else {
      stop("Coordinate df must contain x/y or row/col columns. Found: ", paste(colnames(df), collapse = ", "))
    }
  }

  meta <- obj@meta.data
  meta$cell_id <- obj$cell_id
  safe_pull <- function(vec, keys) vec[keys]

  if (!roi_col %in% colnames(meta)) {
    df[[roi_col]] <- FALSE
  } else {
    roi_vec <- meta[[roi_col]]
    names(roi_vec) <- meta$cell_id
    df[[roi_col]] <- safe_pull(roi_vec, df$cell_id)
  }

  if (!celltype_col %in% colnames(meta)) {
    df[[celltype_col]] <- NA_character_
  } else {
    ct_vec <- meta[[celltype_col]]
    names(ct_vec) <- meta$cell_id
    df[[celltype_col]] <- safe_pull(ct_vec, df$cell_id)
  }

  if (drop_na_celltype) {
    df <- df[!is.na(df[[celltype_col]]), , drop = FALSE]
  }

  if (any(is.na(df[[roi_col]]))) {
    na_frac <- mean(is.na(df[[roi_col]]))
    warning(sprintf(
      "ROI join produced NAs (%.1f%%). This usually means df$cell_id does not match Seurat cell IDs.",
      100 * na_frac
    ))
  }

  computed_xlim <- NULL
  computed_ylim <- NULL
  if (isTRUE(zoom)) {
    if (!is.null(zoom_xlim)) {
      computed_xlim <- zoom_xlim
    }
    if (!is.null(zoom_ylim)) {
      computed_ylim <- zoom_ylim
    }

    if (is.null(computed_xlim) || is.null(computed_ylim)) {
      roi_idx <- which(df[[roi_col]] %in% TRUE)
      if (length(roi_idx) == 0) {
        warning("zoom=TRUE but no ROI cells found (", roi_col, " == TRUE). Zoom ignored.")
      } else {
        x_roi <- df$x[roi_idx]
        y_roi <- df$y[roi_idx]
        xr <- range(x_roi, na.rm = TRUE)
        yr <- range(y_roi, na.rm = TRUE)

        xpad <- diff(xr) * zoom_pad
        ypad <- diff(yr) * zoom_pad
        if (!is.finite(xpad) || xpad == 0) {
          xpad <- 1
        }
        if (!is.finite(ypad) || ypad == 0) {
          ypad <- 1
        }

        if (is.null(computed_xlim)) {
          computed_xlim <- c(xr[1] - xpad, xr[2] + xpad)
        }
        if (is.null(computed_ylim)) {
          computed_ylim <- c(yr[1] - ypad, yr[2] + ypad)
        }
      }
    }
  }

  apply_zoom <- function(p) {
    if (!isTRUE(zoom)) {
      return(p)
    }
    if (is.null(computed_xlim) || is.null(computed_ylim)) {
      return(p)
    }
    p + ggplot2::coord_cartesian(xlim = computed_xlim, ylim = computed_ylim)
  }

  apply_legend_style <- function(p, guide_title = NULL) {
    p +
      ggplot2::theme(
        legend.position = legend_position,
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = ggplot2::element_text(size = legend_title_size),
        legend.text = ggplot2::element_text(size = legend_text_size)
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = guide_title,
          nrow = legend_nrow,
          ncol = legend_ncol,
          override.aes = list(size = legend_point_size, alpha = 1)
        )
      )
  }

  p_roi <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = .data[[roi_col]])) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_manual(values = roi_colors, name = "ROI") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::ggtitle(
      paste0(
        "ROI QC (",
        if (!is.null(obj@misc$xenium_image_name)) obj@misc$xenium_image_name else "image",
        "): ",
        roi_col
      )
    )

  if (isTRUE(grey_outside_roi)) {
    p_celltype <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(
        data = df[df[[roi_col]] == FALSE, ],
        color = "grey80",
        alpha = alpha,
        size = point_size
      ) +
      ggplot2::geom_point(
        data = df[df[[roi_col]] == TRUE, ],
        ggplot2::aes(color = .data[[celltype_col]]),
        alpha = alpha,
        size = point_size
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste0("Cell type: ", celltype_col))
  } else {
    p_celltype <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = .data[[celltype_col]])) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste0("Cell type: ", celltype_col))
  }

  if (!is.null(celltype_colors)) {
    if (is.null(names(celltype_colors))) {
      stop("celltype_colors must be a *named* vector: names = cell type labels, values = colors.")
    }

    present <- unique(df[[celltype_col]])
    present <- present[!is.na(present)]

    missing <- setdiff(present, names(celltype_colors))
    if (length(missing) > 0) {
      msg <- paste0("Missing colors for: ", paste(missing, collapse = ", "))
      if (isTRUE(strict_celltype_colors)) {
        stop(msg)
      } else {
        warning(msg)
      }
    }

    p_celltype <- p_celltype +
      ggplot2::scale_color_manual(values = celltype_colors, name = celltype_col, drop = FALSE)
  }

  if (isTRUE(flip_y)) {
    p_roi <- p_roi + ggplot2::scale_y_reverse()
    p_celltype <- p_celltype + ggplot2::scale_y_reverse()
  }

  if (isTRUE(zoom)) {
    if (isTRUE(zoom_fixed)) {
      p_roi <- apply_zoom(p_roi)
      p_celltype <- apply_zoom(p_celltype)
    } else {
      p_roi <- apply_zoom(p_roi)
    }
  }

  p_roi <- apply_legend_style(p_roi, guide_title = "ROI")
  p_celltype <- apply_legend_style(p_celltype, guide_title = celltype_col)

  if (show == "roi") {
    return(list(df = df, p = p_roi, xlim = computed_xlim, ylim = computed_ylim))
  }
  if (show == "celltype") {
    return(list(df = df, p = p_celltype, xlim = computed_xlim, ylim = computed_ylim))
  }

  combined <- patchwork::wrap_plots(p_roi, p_celltype, nrow = 1)
  if (isTRUE(collect_legends)) {
    combined <- combined +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = legend_position)
  }

  list(df = df, p_roi = p_roi, p_celltype = p_celltype, p = combined,
       xlim = computed_xlim, ylim = computed_ylim)
}

#' Plot Cell Type Distribution Inside ROI
#'
#' Creates a horizontal bar plot of cell type counts or proportions among ROI cells.
#'
#' @param df Data frame with ROI and cell type columns
#' @param roi_col ROI indicator column name
#' @param celltype_col Cell type column name
#' @param celltype_colors Optional named cell type color vector
#' @param metric One of "count" or "proportion"
#' @param order_by One of "value", "palette", "none"
#' @param legend Logical; show legend
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Optional y-axis label
#' @param labels One of "none", "count", "percent", "both"
#' @param label_size Label text size
#' @param label_hjust Label horizontal adjustment
#' @param label_accuracy Percent label rounding accuracy
#' @param label_pad_mult Extra axis expansion for labels
#' @return ggplot object
#' @export
plot_roi_celltype_dist <- function(
    df,
    roi_col,
    celltype_col,
    celltype_colors = NULL,
    metric = c("count", "proportion"),
    order_by = c("value", "palette", "none"),
    legend = FALSE,
    title = "Cell type distribution inside ROI",
    xlab = "Cell type",
    ylab = NULL,
    labels = c("none", "count", "percent", "both"),
    label_size = 3,
    label_hjust = -0.1,
    label_accuracy = 0.1,
    label_pad_mult = 0.15
) {
  # Auto-apply brca atlas color palette
  if (is.null(celltype_colors) && celltype_col == "singleR.predicted.id.brcaAtlas") {
    celltype_colors <- brca_atlas_colors()
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install dplyr: install.packages('dplyr')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Please install scales: install.packages('scales')")
  }

  metric <- match.arg(metric)
  order_by <- match.arg(order_by)
  labels <- match.arg(labels)

  if (!roi_col %in% colnames(df)) {
    stop("roi_col not found in df: ", roi_col)
  }
  if (!celltype_col %in% colnames(df)) {
    stop("celltype_col not found in df: ", celltype_col)
  }

  d <- df |>
    dplyr::filter(.data[[roi_col]] == TRUE) |>
    dplyr::count(.data[[celltype_col]], name = "n")

  if (nrow(d) == 0) {
    stop("No rows remain after filtering ", roi_col, " == TRUE. Check ROI column or data.")
  }

  d <- d |>
    dplyr::mutate(
      prop = n / sum(n),
      value = if (metric == "proportion") prop else n
    )

  if (is.null(ylab)) {
    ylab <- if (metric == "proportion") "Proportion of ROI cells" else "Number of ROI cells"
  }

  if (order_by == "value") {
    d <- d |>
      dplyr::mutate(.ct = stats::reorder(.data[[celltype_col]], value))
  } else if (order_by == "palette") {
    if (is.null(celltype_colors) || is.null(names(celltype_colors))) {
      stop("order_by = 'palette' requires celltype_colors as a *named* vector.")
    }
    d <- d |>
      dplyr::mutate(.ct = factor(.data[[celltype_col]], levels = names(celltype_colors)))
  } else {
    d <- d |>
      dplyr::mutate(.ct = .data[[celltype_col]])
  }

  if (labels == "count") {
    d$lbl <- as.character(d$n)
  } else if (labels == "percent") {
    d$lbl <- scales::percent(d$prop, accuracy = label_accuracy)
  } else if (labels == "both") {
    d$lbl <- paste0(d$n, " (", scales::percent(d$prop, accuracy = label_accuracy), ")")
  } else {
    d$lbl <- NA_character_
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .ct, y = value, fill = .ct)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title
    )

  if (!is.null(celltype_colors)) {
    if (is.null(names(celltype_colors))) {
      stop("celltype_colors must be a *named* vector: names = cell type labels, values = colors.")
    }
    p <- p + ggplot2::scale_fill_manual(values = celltype_colors, name = "Cell type", drop = FALSE)
  } else {
    p <- p + ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell type"))
  }

  if (metric == "proportion") {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent,
      expand = ggplot2::expansion(mult = c(0, label_pad_mult))
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, label_pad_mult))
    )
  }

  if (labels != "none") {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = lbl),
        hjust = label_hjust,
        size = label_size
      )
  }

  if (!isTRUE(legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}

#' Plot Stacked Cell Type Composition Inside vs Outside ROI
#'
#' Creates stacked bars comparing cell type composition for inside and outside ROI.
#'
#' @param df Data frame with ROI and cell type columns
#' @param roi_col ROI indicator column name
#' @param celltype_col Cell type column name
#' @param celltype_colors Optional named cell type color vector
#' @param metric One of "count" or "proportion"
#' @param order_by One of "inside", "total", "palette", "none"
#' @param inside_label Label for inside ROI group
#' @param outside_label Label for outside ROI group
#' @param legend Logical; show legend
#' @param legend_title Legend title
#' @param legend_position Legend position
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Optional y-axis label
#' @return ggplot object
#' @export
plot_roi_inout_celltype_stack <- function(
    df,
    roi_col,
    celltype_col,
    celltype_colors = NULL,
    metric = c("count", "proportion"),
    order_by = c("inside", "total", "palette", "none"),
    inside_label = "Inside ROI",
    outside_label = "Outside ROI",
    legend = TRUE,
    legend_title = "Cell type",
    legend_position = "bottom",
    title = "Cell type distribution: Inside vs Outside ROI",
    xlab = NULL,
    ylab = NULL
) {
  # Auto-apply brca atlas color palette
  if (is.null(celltype_colors) && celltype_col == "singleR.predicted.id.brcaAtlas") {
    celltype_colors <- brca_atlas_colors()
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install dplyr: install.packages('dplyr')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Please install scales: install.packages('scales')")
  }

  metric <- match.arg(metric)
  order_by <- match.arg(order_by)

  if (!roi_col %in% colnames(df)) {
    stop("roi_col not found in df: ", roi_col)
  }
  if (!celltype_col %in% colnames(df)) {
    stop("celltype_col not found in df: ", celltype_col)
  }

  d <- df |>
    dplyr::mutate(
      roi_status = ifelse(.data[[roi_col]] %in% TRUE, inside_label, outside_label)
    ) |>
    dplyr::count(roi_status, .data[[celltype_col]], name = "n")

  # Factor cell types turn into integer codes when combined with c() below
  d[[celltype_col]] <- as.character(d[[celltype_col]])
  d$roi_status <- factor(d$roi_status, levels = c(inside_label, outside_label))

  d <- d |>
    dplyr::group_by(roi_status) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  if (metric == "proportion") {
    d$value <- d$prop
    if (is.null(ylab)) {
      ylab <- "Proportion of cells"
    }
  } else {
    d$value <- d$n
    if (is.null(ylab)) {
      ylab <- "Number of cells"
    }
  }
  if (is.null(xlab)) {
    xlab <- NULL
  }

  ct_vals <- d[[celltype_col]]
  if (order_by == "inside") {
    ct_order <- d |>
      dplyr::filter(roi_status == inside_label) |>
      dplyr::arrange(dplyr::desc(n)) |>
      dplyr::pull(.data[[celltype_col]])

    remaining <- setdiff(unique(ct_vals), ct_order)
    ct_order <- c(ct_order, remaining)
    d[[celltype_col]] <- factor(ct_vals, levels = ct_order)
  } else if (order_by == "total") {
    ct_order <- d |>
      dplyr::group_by(.data[[celltype_col]]) |>
      dplyr::summarise(total_n = sum(n), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(total_n)) |>
      dplyr::pull(.data[[celltype_col]])

    d[[celltype_col]] <- factor(ct_vals, levels = ct_order)
  } else if (order_by == "palette") {
    if (is.null(celltype_colors) || is.null(names(celltype_colors))) {
      stop("order_by = 'palette' requires celltype_colors as a *named* vector.")
    }
    d[[celltype_col]] <- factor(ct_vals, levels = names(celltype_colors))
  } else {
    d[[celltype_col]] <- ct_vals
  }

  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = roi_status, y = value, fill = .data[[celltype_col]])
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title,
      fill = legend_title
    )

  if (metric == "proportion") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent)
  }

  if (!is.null(celltype_colors)) {
    if (is.null(names(celltype_colors))) {
      stop("celltype_colors must be a *named* vector: names = cell type labels, values = colors.")
    }
    p <- p + ggplot2::scale_fill_manual(values = celltype_colors, drop = FALSE)
  }

  if (!isTRUE(legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = legend_position)
  }

  return(p)
}

#' Plot Dodged Cell Type Comparison Inside vs Outside ROI
#'
#' Creates a horizontal dodged bar plot comparing inside and outside ROI counts
#' or proportions by cell type.
#'
#' @param df Data frame with ROI and cell type columns
#' @param roi_col ROI indicator column name
#' @param celltype_col Cell type column name
#' @param metric One of "count" or "proportion"
#' @param roi_colors Named colors for `FALSE` and `TRUE`
#' @param inside_label Label for inside ROI group
#' @param outside_label Label for outside ROI group
#' @param order_by One of "total", "inside", "none"
#' @param top_n Optional number of top cell types to keep
#' @param top_n_by One of "total" or "inside"
#' @param other_label Label for collapsed non-top groups
#' @param keep_other Logical; keep collapsed other group
#' @param legend Logical; show legend
#' @param legend_title Legend title
#' @param legend_position Legend position
#' @param title Optional plot title
#' @param xlab X-axis label
#' @param ylab Optional y-axis label
#' @param labels One of "none", "count", "percent", "both"
#' @param label_size Label text size
#' @param label_hjust Label horizontal adjustment
#' @param label_pad_mult Extra axis expansion for labels
#' @param label_accuracy Percent label rounding accuracy
#' @return ggplot object
#' @export
plot_roi_inout_celltype_dodge <- function(
    df,
    roi_col,
    celltype_col,
    metric = c("count", "proportion"),
    roi_colors = c(`FALSE` = "grey80", `TRUE` = "red"),
    inside_label = "Inside ROI",
    outside_label = "Outside ROI",
    order_by = c("total", "inside", "none"),
    top_n = NULL,
    top_n_by = c("total", "inside"),
    other_label = "Other",
    keep_other = TRUE,
    legend = TRUE,
    legend_title = "ROI status",
    legend_position = "right",
    title = NULL,
    xlab = "Cell type",
    ylab = NULL,
    labels = c("none", "count", "percent", "both"),
    label_size = 3,
    label_hjust = -0.1,
    label_pad_mult = 0.15,
    label_accuracy = 0.1
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install dplyr")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Please install scales")
  }

  metric <- match.arg(metric)
  order_by <- match.arg(order_by)
  labels <- match.arg(labels)
  top_n_by <- match.arg(top_n_by)

  if (!roi_col %in% colnames(df)) {
    stop("roi_col not found in df: ", roi_col)
  }
  if (!celltype_col %in% colnames(df)) {
    stop("celltype_col not found in df: ", celltype_col)
  }
  if (is.null(names(roi_colors)) || !all(c("FALSE", "TRUE") %in% names(roi_colors))) {
    stop("roi_colors must be named with 'FALSE' and 'TRUE' (e.g. c(`FALSE`='grey80', `TRUE`='red')).")
  }
  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1)) {
    stop("top_n must be NULL or a single integer >= 1.")
  }

  d <- df |>
    dplyr::mutate(
      roi_status = ifelse(.data[[roi_col]] %in% TRUE, inside_label, outside_label)
    ) |>
    dplyr::count(roi_status, .data[[celltype_col]], name = "n")

  # Factor cell types turn into integer codes when combined with c() below
  d[[celltype_col]] <- as.character(d[[celltype_col]])
  d$roi_status <- factor(d$roi_status, levels = c(inside_label, outside_label))

  if (!is.null(top_n)) {
    rank_tbl <- d |>
      dplyr::group_by(.data[[celltype_col]]) |>
      dplyr::summarise(
        total_n = sum(n),
        inside_n = sum(n[roi_status == inside_label]),
        .groups = "drop"
      )

    top_types <- if (top_n_by == "inside") {
      rank_tbl |>
        dplyr::arrange(dplyr::desc(inside_n), dplyr::desc(total_n)) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::pull(.data[[celltype_col]])
    } else {
      rank_tbl |>
        dplyr::arrange(dplyr::desc(total_n)) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::pull(.data[[celltype_col]])
    }

    if (isTRUE(keep_other)) {
      d <- d |>
        dplyr::mutate(
          ct_grp = dplyr::if_else(.data[[celltype_col]] %in% top_types,
                                  as.character(.data[[celltype_col]]),
                                  other_label)
        ) |>
        dplyr::group_by(roi_status, ct_grp) |>
        dplyr::summarise(n = sum(n), .groups = "drop")

      names(d)[names(d) == "ct_grp"] <- celltype_col

      top_types <- c(top_types, other_label)
    } else {
      d <- d |>
        dplyr::filter(.data[[celltype_col]] %in% top_types)
    }
  }

  d <- d |>
    dplyr::group_by(roi_status) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  if (metric == "proportion") {
    d$value <- d$prop
    if (is.null(ylab)) {
      ylab <- "Proportion of cells"
    }
    if (is.null(title)) {
      title <- "Cell type proportions by cell type inside vs outside ROI"
    }
  } else {
    d$value <- d$n
    if (is.null(ylab)) {
      ylab <- "Number of cells"
    }
    if (is.null(title)) {
      title <- "Cell counts by cell type inside vs outside ROI"
    }
  }

  if (order_by == "total") {
    ct_order <- d |>
      dplyr::group_by(.data[[celltype_col]]) |>
      dplyr::summarise(total_n = sum(n), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(total_n)) |>
      dplyr::pull(.data[[celltype_col]])
    d[[celltype_col]] <- factor(d[[celltype_col]], levels = ct_order)
  } else if (order_by == "inside") {
    ct_order <- d |>
      dplyr::filter(roi_status == inside_label) |>
      dplyr::arrange(dplyr::desc(n)) |>
      dplyr::pull(.data[[celltype_col]])
    remaining <- setdiff(unique(as.character(d[[celltype_col]])), ct_order)
    d[[celltype_col]] <- factor(as.character(d[[celltype_col]]), levels = c(ct_order, remaining))
  }

  if (labels == "count") {
    d$lbl <- as.character(d$n)
  } else if (labels == "percent") {
    d$lbl <- scales::percent(d$prop, accuracy = label_accuracy)
  } else if (labels == "both") {
    d$lbl <- paste0(d$n, " (", scales::percent(d$prop, accuracy = label_accuracy), ")")
  } else {
    d$lbl <- NA_character_
  }

  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data[[celltype_col]], y = value, fill = roi_status)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = xlab, y = ylab, title = title)

  if (metric == "proportion") {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent,
      expand = ggplot2::expansion(mult = c(0, label_pad_mult))
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, label_pad_mult))
    )
  }

  p <- p + ggplot2::scale_fill_manual(
    values = stats::setNames(
      c(roi_colors[["TRUE"]], roi_colors[["FALSE"]]),
      c(inside_label, outside_label)
    ),
    name = legend_title
  )

  if (labels != "none") {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = lbl),
        position = ggplot2::position_dodge(width = 0.8),
        hjust = label_hjust,
        size = label_size
      )
  }

  if (!isTRUE(legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = legend_position)
  }

  p
}
