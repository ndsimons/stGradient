# Plot DBSCAN-based cancer region polygons

Uses DBSCAN clustering on cancer cell coordinates to identify spatial
regions, creates concave hull polygons around each cluster, and
optionally visualizes them.

## Usage

``` r
plot_dbscan_polygons(
  seurat_obj,
  fov = "fov",
  eps = 35,
  minPts = 20,
  concavity = 2,
  length_threshold = 0,
  sample_name = NULL,
  plot = TRUE,
  ncol = 2
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates

- fov:

  Name of the field of view (default: "fov")

- eps:

  DBSCAN epsilon parameter - max distance between points (default: 35)

- minPts:

  DBSCAN minimum points per cluster (default: 20)

- concavity:

  Concave hull concavity parameter, lower = tighter fit (default: 2)

- length_threshold:

  Edge length threshold for concave hull (default: 0)

- sample_name:

  Optional sample name prefixed to the plot titles

- plot:

  Whether to generate visualization plots (default: TRUE)

- ncol:

  Number of plot columns: 2 places the plots side by side, 1 stacks them
  vertically (default: 2)

## Value

sf object containing polygon geometries for each cancer region

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
polygons <- plotDBSCANpolygons(seurat_obj, sample_name = "Patient1")

# Custom parameters without plotting
polygons <- plotDBSCANpolygons(
  seurat_obj, 
  eps = 50, 
  minPts = 30,
  sample_name = "Patient1",
  plot = FALSE
)
} # }
```
