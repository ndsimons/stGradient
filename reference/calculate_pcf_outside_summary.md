# Calculate Outside-Polygon Pair Correlation Summary

Computes pair correlation function (PCF) values for all requested
cell-type pairs using only cells outside the tumor polygons.

## Usage

``` r
calculate_pcf_outside_summary(
  seurat_obj,
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
  bw = 10
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- polygons:

  sf polygon object representing ROIs.

- distance:

  Target distance (microns) at which to extract the PCF value. g is
  evaluated at exactly this distance.

- min_cells:

  Minimum number of cells required per cell type.

- celltype_col:

  Metadata column containing cell type labels.

- fov_name:

  Field of view/image name passed to \`Seurat::GetTissueCoordinates\`.

- include_self_pairs:

  Logical; include same-type pairs.

- ordered_pairs:

  Logical; if \`TRUE\`, computes ordered pairs (A,B and B,A). If
  \`FALSE\`, computes each pair once.

- window_obj:

  Optional \`spatstat.geom::owin\` tissue window; tumor polygons are
  subtracted from it.

- show_progress:

  Logical; show a text progress bar.

- tissue_window:

  Tissue outline used when \`window_obj\` is \`NULL\`: \`"concave"\`
  (default; concave hull of all cells), \`"convex"\`, or \`"bbox"\`.
  Tumor polygons are always removed from the window.

- bw:

  Kernel bandwidth (microns) passed to \`pcfcross\`. Use the same value
  as in \[plot_pcf_outside_pair_envelopes()\] so heatmap values match
  the envelope curves. \`NULL\` uses the \`spatstat\` default.

## Value

Data frame with columns \`type1\`, \`type2\`, \`g_value\`, \`n_type1\`,
\`n_type2\`.
