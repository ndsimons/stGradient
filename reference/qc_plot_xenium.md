# QC Visualization for Xenium Spatial Data

Generates ROI and cell type overlays from spatial coordinates and
metadata.

## Usage

``` r
qc_plot_xenium(
  obj,
  roi_col = "inside_polygon",
  celltype_col,
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
)
```

## Arguments

- obj:

  Seurat object with spatial coordinates

- roi_col:

  Metadata column indicating ROI membership

- celltype_col:

  Metadata column containing cell type annotation

- celltype_colors:

  Optional named vector mapping cell types to colors

- strict_celltype_colors:

  Logical; error if present cell types are missing colors

- image:

  Optional image/FOV name

- coords_df:

  Optional pre-built coordinate data.frame with \`cell_id\`

- point_size:

  Point size for scatter layers

- alpha:

  Point alpha

- flip_y:

  Logical; reverse y-axis

- show:

  One of "both", "roi", "celltype"

- roi_colors:

  Named vector of colors for \`FALSE\` and \`TRUE\` ROI values

- drop_na_celltype:

  Logical; drop rows with missing cell type

- grey_outside_roi:

  Logical; draw outside ROI cells in gray on celltype panel

- zoom:

  Logical; zoom to ROI bounding box or manual limits

- zoom_pad:

  Fractional padding around ROI bounding box

- zoom_fixed:

  Logical; apply same zoom limits to both panels

- zoom_xlim:

  Optional numeric vector \`c(xmin, xmax)\`

- zoom_ylim:

  Optional numeric vector \`c(ymin, ymax)\`

- legend_position:

  Legend position passed to ggplot2 theme

- legend_point_size:

  Point size in legend keys

- legend_title_size:

  Legend title text size

- legend_text_size:

  Legend text size

- legend_nrow:

  Optional legend rows

- legend_ncol:

  Optional legend columns

- collect_legends:

  Logical; collect legends when \`show = "both"\`

## Value

A list containing processed data and one or more ggplot objects
