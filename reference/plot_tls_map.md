# Plot Called TLS on the Tissue

Spatial map of all cells (grey) with cells in called TLS colored by TLS
ID, B cells within TLS in red, and TLS IDs labeled at their centroids.

## Usage

``` r
plot_tls_map(
  seurat_obj,
  fov_name = "fov",
  celltype_col = "first_type",
  b_types = c("B_Cells"),
  point_size = 0.5
)
```

## Arguments

- seurat_obj:

  Seurat object processed with \[call_tls()\].

- fov_name:

  Field of view/image name (default: "fov").

- celltype_col:

  Metadata column with cell type labels.

- b_types:

  Cell type labels to highlight as B cells.

- point_size:

  Point size for TLS cells.

## Value

ggplot object.
