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
  celltype_col = "singleR.predicted.id.brcaAtlas",
  cancer_labels = c("Cancer LumA SC", "Cancer LumB SC", "Cancer Her2 SC")
)
```

See: [Loading & Setup](loading-and-setup.md)

## Step 2: Define Tumor Regions

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
  group_by = "cellType",
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

## Advanced Analyses

Once the basic distance profile is computed, you can extend into:

- **PCF Analysis** — characterize cell-type spatial co-occurrence
  outside the tumor ([PCF Analysis](pcf-analysis.md))
- **Gene-Distance Modeling** — identify genes whose expression varies
  with boundary distance ([Gene-Distance
  Analysis](gene-distance-analysis.md))
