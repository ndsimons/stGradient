# stGradient Workflow

## Overview

This vignette walks through a complete stGradient analysis from loading
data to final visualization. For detailed documentation on each step,
see the topic-specific vignettes linked below.

## Step 1: Load Data

``` r

library(stGradient)

seurat_obj <- load_xenium_seurat(
  rds_path = "path/to/sample.rds",
  celltype_col = "first_type",
  cancer_labels = c("DCIS 1", "DCIS 2", "Invasive_Tumor", "Prolif_Invasive_Tumor")
)
```

See: [Loading & Setup](loading-and-setup.md)

## Step 2: Define ROIs

Option A — Interactive (Shiny app):

``` r

launch_stGradient()
seurat_obj <- add_polygon_membership(seurat_obj, drawn_polygons)
```

Option B — Automated (DBSCAN):

``` r

seurat_obj <- add_dbscan_polygons(seurat_obj)
auto_polygons <- get_cancer_regions(seurat_obj)
```

See: [Region Selection](region-selection.md)

## Step 3: Calculate Boundary Distances

``` r

seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, auto_polygons)
```

See: [Distance Metrics](distance-metrics.md)

## Step 4: Distance Profile Analysis

``` r

results <- calculate_distance_profile(
  seurat_obj = seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "first_type",
  inside_mode = "all"
)
```

## Step 5: Find Critical Distances

``` r

critical <- find_critical_distance(results$profile, "transcriptomic_distance")
cat("Maximum separation at:", critical$max_distance, "μm\n")
```

## Step 6: Visualize

``` r

plot_all_metrics(results$profile)
plot_celltype_heatmap(results$celltype_contributions)
plot_celltype_proportions(results$celltype_contributions, top_n = 5)
```

See: [Visualization](visualization.md)

## Step 7: PCF Analysis

Characterize cell-type spatial co-occurrence outside the tumor with pair
correlation functions. Requires the `spatstat.geom` and
`spatstat.explore` packages.

``` r

pcf_summary <- calculate_pcf_outside_summary(
  seurat_obj = seurat_obj,
  polygons = auto_polygons,
  distance = 200,
  min_cells = 10,
  celltype_col = "first_type",
  include_self_pairs = TRUE,
  ordered_pairs = FALSE
)
plot_pcf_outside_heatmap(pcf_summary, midpoint = 1)

# Monte Carlo envelopes for every pair, saved as one PNG per pair
env_list <- plot_pcf_outside_pair_envelopes(
  seurat_obj = seurat_obj,
  polygons = auto_polygons,
  celltype_col = "first_type",
  min_cells = 10,
  nsim = 39,
  file = "pcf_outside_all_pairs.png"
)
plot_pcf_pair_gg(env_list, "TypeA", "TypeB")
```

g(r) \> 1 indicates clustering and g(r) \< 1 indicates avoidance.

See: [PCF Analysis](pcf-analysis.md)

## Step 8: Gene-Distance Analysis

Identify genes whose expression varies with signed distance from the
tumor boundary.

``` r

seurat_obj <- compute_distance_from_boundary(seurat_obj, fov_name = "fov")

dist_res <- find_distance_associated_genes(
  seurat_obj,
  assay = "Xenium",
  focus = "peritumoral",
  min_cells_expressing = 10,
  n_distance_bins = 10,
  max_distance = 200,
  p_adj_threshold = 0.05
)

dist_plots <- plot_distance_genes(dist_res, top_n = 12, sample_label = "Sample_01")
dist_plots$volcano
dist_plots$trends
```

See: [Gene-Distance Analysis](gene-distance-analysis.md)

## Step 9: TLS Analysis

Call tertiary lymphoid structures from B and T cell aggregates, then
validate them with a TLS gene signature score. Set the cell type column
and labels to match your annotation.

``` r

seurat_obj <- call_tls(
  seurat_obj,
  fov_name = "fov",
  celltype_col = "first_type",
  b_types = c("B_Cells"),
  t_types = c("CD4+_T_Cells", "CD8+_T_Cells"),
  dc_types = c("LAMP3+_DCs", "IRF7+_DCs"),
  eps = 50,
  minPts = 20
)
get_tls_summary(seurat_obj)
plot_tls_map(seurat_obj, celltype_col = "first_type")
plot_tls_composition(seurat_obj, celltype_col = "first_type")

seurat_obj <- add_tls_score(seurat_obj, assay = "Xenium")
plot_tls_score(seurat_obj, highlight = "outline")

tls_stats <- summarize_tls_score(seurat_obj)
tls_stats$by_tls
plot_tls_score_comparison(seurat_obj, type = "violin")
```

See: [TLS Analysis](tls-analysis.md)
