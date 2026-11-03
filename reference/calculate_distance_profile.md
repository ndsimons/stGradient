# Calculate Distance Profile

Calculates genomic separation metrics at varying distances from boundary

## Usage

``` r
calculate_distance_profile(
  seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by,
  inside_mode = "all"
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial data

- distance_thresholds:

  Vector of distances to test (in microns)

- metric:

  Which metrics to calculate ("all", "transcriptomic", "composition",
  "de_strength")

- group_by:

  Metadata column for cell type composition

- inside_mode:

  "all" = all cells inside polygon, "distance" = only cells within
  threshold

## Value

List with profile data frame and celltype contributions
