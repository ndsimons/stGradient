# Compute Signed Distance From Tumor Boundary

Computes per-cell signed distance to the nearest cancer region polygon,
using the same per-region polygons and definition as
\[calculate_boundary_distances()\]. \`dist_to_boundary\` is therefore
identical to \`boundary_distance\` computed from the same polygons, and
stroma between separate tumor nests is classed as outside the tumor.

## Usage

``` r
compute_distance_from_boundary(
  seurat_obj,
  polygons = NULL,
  fov_name = "fov",
  polygon_col = NULL,
  intratumoral_label = NULL,
  concavity = NULL
)
```

## Arguments

- seurat_obj:

  Seurat object.

- polygons:

  sf polygon object. Defaults to the polygons stored by
  \[add_dbscan_polygons()\], retrieved with \[get_cancer_regions()\].

- fov_name:

  Field of view/image name for coordinate extraction.

- polygon_col, intratumoral_label, concavity:

  Deprecated and ignored. The boundary used to be a single concave hull
  around all intratumoral cells, which merged separate tumor nests and
  labelled the stroma between them as intratumoral.

## Value

Seurat object with \`dist_to_boundary\` and \`abs_dist_to_boundary\`
metadata.
