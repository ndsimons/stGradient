# Optimize DBSCAN Parameters

Performs a grid search over eps and minPts for DBSCAN clustering of
cancer cell coordinates. Each combination is scored by average
silhouette score and by the fraction of cancer cells labeled noise.
Combinations with noise below \`noise_pct\` are turned into concave hull
polygons and scored by coverage (fraction of cancer cells inside a
polygon) and purity (fraction of cells inside polygons that are cancer),
which are combined into an F-beta score.

## Usage

``` r
optimize_dbscan_params(
  seurat_obj,
  noise_pct = 0.15,
  b = 0.5,
  plot = TRUE,
  fov_name = "fov",
  cancer_col = "is_cancer",
  eps_range = seq(10, 100, by = 5),
  minPts_range = c(10, 15, 20, 30, 50),
  concavity = 2,
  length_threshold = 0,
  max_sample = 46000,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates

- noise_pct:

  Maximum fraction of cancer cells labeled noise for a combination to be
  scored by coverage and purity (default: 0.15)

- b:

  Beta for the F-beta score. Values below 1 weight purity more than
  coverage; values above 1 weight coverage more (default: 0.5)

- plot:

  Print heatmaps of coverage, purity and F-beta for the scored
  combinations, with the best combination outlined in red on the F-beta
  heatmap (default: TRUE)

- fov_name:

  Name of the field of view (default: "fov")

- cancer_col:

  Metadata column indicating cancer cells (default: "is_cancer")

- eps_range:

  Numeric vector of eps values to test (default: seq(10, 100, by = 5))

- minPts_range:

  Numeric vector of minPts values to test (default: c(10, 15, 20, 30,
  50))

- concavity:

  Concave hull concavity parameter used to build polygons (default: 2)

- length_threshold:

  Edge length threshold for concave hull (default: 0)

- max_sample:

  Maximum number of non-noise points to sample for silhouette
  calculation to control memory usage (default: 46000)

- verbose:

  Print progress messages (default: TRUE)

## Value

A data frame with columns eps, minPts, sil_score, n_clusters, noise_pct,
n_regions, cancer_covered, purity and fbeta, sorted by descending fbeta
(then silhouette score). The best combination is in the first row.
Combinations above \`noise_pct\` have \`NA\` coverage, purity and fbeta.
When \`plot = TRUE\`, the heatmap is also stored in \`attr(result,
"plot")\`.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- optimize_dbscan_params(seurat_obj, noise_pct = 0.15, b = 0.5, plot = TRUE)
head(results)

# Use best parameters
best <- results[1, ]
seurat_obj <- add_dbscan_polygons(seurat_obj, eps = best$eps, minPts = best$minPts)
} # }
```
