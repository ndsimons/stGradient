# Automatically detect cancer cell regions and create polygons

Uses DBSCAN clustering on cancer cell coordinates to identify spatial
regions, then creates concave hull polygons around each cluster. Stores
the polygons in \`seurat_obj@misc\$cancer_regions\[\[fov_name\]\]\`
(retrieve them with \[get_cancer_regions()\]) and adds an
\`inside_polygon\` metadata column.

## Usage

``` r
add_dbscan_polygons(
  seurat_obj,
  fov_name = "fov",
  cancer_col = "is_cancer",
  eps = 35,
  minPts = 20,
  concavity = 2,
  length_threshold = 0
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates

- fov_name:

  Name of the field of view (default: "fov")

- cancer_col:

  Metadata column indicating cancer cells (default: "is_cancer")

- eps:

  DBSCAN epsilon parameter - max distance between points (default: 35)

- minPts:

  DBSCAN minimum points per cluster (default: 20)

- concavity:

  Concave hull concavity parameter, lower = tighter fit (default: 2)

- length_threshold:

  Edge length threshold for concave hull (default: 0)

## Value

Seurat object with cancer region polygons stored in \`@misc\` and
\`inside_polygon\` metadata

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- add_dbscan_polygons(seurat_obj, eps = 50, minPts = 30)
} # }
```
