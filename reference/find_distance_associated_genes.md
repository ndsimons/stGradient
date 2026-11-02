# Find Distance-Associated Genes

Fits per-gene linear models and Spearman correlations against distance
from the tumor boundary.

## Usage

``` r
find_distance_associated_genes(
  seurat_obj,
  assay = "Xenium",
  min_cells_expressing = 10,
  n_distance_bins = 10,
  focus = c("peritumoral", "intratumoral", "both"),
  max_dist_quantile = 0.95,
  max_distance = 200,
  p_adj_threshold = 0.05
)
```

## Arguments

- seurat_obj:

  Seurat object with \`dist_to_boundary\` metadata.

- assay:

  Assay name. If its \`data\` layer is empty,
  \`Seurat::NormalizeData()\` is run with default settings on the
  analyzed cells.

- min_cells_expressing:

  Minimum expressing cells for a gene to be tested.

- n_distance_bins:

  Number of bins for trend summaries.

- focus:

  Region to analyze: \`peritumoral\`, \`intratumoral\`, or \`both\`.

- max_dist_quantile:

  Distance quantile used to trim extreme outliers.

- max_distance:

  Maximum distance from the boundary to include (for \`intratumoral\`,
  maximum absolute distance). Cells farther away are excluded from model
  fitting and binning. Combined with \`max_dist_quantile\`: the smaller
  of the two cutoffs is used. Set to \`NULL\` to rely on the quantile
  alone.

- p_adj_threshold:

  BH adjusted p-value threshold for significant genes.

## Value

List containing model results, binned expression, and metadata.
