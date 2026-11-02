# Find Genes Associated With Distance From T Cell Hotspots

Find Genes Associated With Distance From T Cell Hotspots

## Usage

``` r
find_hotspot_distance_genes(
  seurat_obj,
  hotspot_result,
  assay = "Xenium",
  focus = c("perihotspot", "intrahotspot", "both"),
  min_cells_expressing = 10,
  n_distance_bins = 10,
  max_dist_quantile = 0.95,
  p_adj_threshold = 0.05,
  exclude_tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
  celltype_col = "singleR.predicted.id.brcaAtlas"
)
```

## Arguments

- seurat_obj:

  Seurat object.

- hotspot_result:

  Output from \`define_tcell_hotspots\`.

- assay:

  Assay name. If its \`data\` layer is empty,
  \`Seurat::NormalizeData()\` is run with default settings on the
  analyzed cells.

- focus:

  Region to analyze relative to hotspot: \`perihotspot\`,
  \`intrahotspot\`, or \`both\`.

- min_cells_expressing:

  Minimum expressing cells for a gene to be tested.

- n_distance_bins:

  Number of bins for trend summaries.

- max_dist_quantile:

  Distance quantile used to trim extreme outliers.

- p_adj_threshold:

  BH adjusted p-value threshold.

- exclude_tcell_labels:

  Optional labels to exclude from analysis.

- celltype_col:

  Metadata cell type column.

## Value

List of model results and binned trend data.
