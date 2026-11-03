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
  celltype_col = "singleR.predicted.id.brcaAtlas",
  fov_name = "fov",
  include_self_pairs = TRUE,
  ordered_pairs = FALSE,
  window_obj = NULL,
  show_progress = interactive()
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- polygons:

  sf polygon object representing ROIs.

- distance:

  Target distance at which to extract PCF value.

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

  Optional \`spatstat.geom::owin\` window.

- show_progress:

  Logical; show a text progress bar.

## Value

Data frame with columns \`type1\`, \`type2\`, \`g_value\`, \`n_type1\`,
\`n_type2\`.
