# Plot Stacked Cell Type Composition Inside vs Outside ROI

Creates stacked bars comparing cell type composition for inside and
outside ROI.

## Usage

``` r
plot_roi_inout_celltype_stack(
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
)
```

## Arguments

- df:

  Data frame with ROI and cell type columns

- roi_col:

  ROI indicator column name

- celltype_col:

  Cell type column name

- celltype_colors:

  Optional named cell type color vector. If \`NULL\`, a built-in
  qualitative palette is assigned to the cell types present.

- metric:

  One of "count" or "proportion"

- order_by:

  One of "inside", "total", "palette", "none"

- inside_label:

  Label for inside ROI group

- outside_label:

  Label for outside ROI group

- legend:

  Logical; show legend

- legend_title:

  Legend title

- legend_position:

  Legend position

- title:

  Plot title

- xlab:

  X-axis label

- ylab:

  Optional y-axis label

## Value

ggplot object
