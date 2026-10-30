# stGradient ![](reference/figures/stGradientLogo.png)

### R toolkit for distance-based analyses of spatial transcriptomics data

stGradient provides tools to interactively select tumor regions and
analyze genomic separation at varying distances from tumor boundaries.
Designed for spatial transcriptomics analysis, particularly for
understanding tumor microenvironment architecture.

## Installation

``` r

devtools::install_github("ndsimons/stGradient")
```

## Getting Started

A [**vignette**](articles/stGradient-workflow.md) is available that
covers the core workflow from loading data through visualization.

Additional documentation organized by topic:

| Vignette | Description |
|----|----|
| [Loading & Setup](articles/loading-and-setup.md) | Load Xenium data, launch the interactive Shiny selector |
| [Region Selection](articles/region-selection.md) | Interactive polygon drawing and DBSCAN-based auto-detection |
| [Distance Metrics](articles/distance-metrics.md) | Boundary distances, transcriptomic distance, composition, DE strength |
| [PCF Analysis](articles/pcf-analysis.md) | Pair correlation functions outside tumor polygons |
| [Gene-Distance Analysis](articles/gene-distance-analysis.md) | Model gene expression as a function of distance from boundaries |
| [Visualization](articles/visualization.md) | QC plots, ROI composition, distance profile plots |

## Quick Start

``` r

library(stGradient)

# Load a Xenium Seurat object
seurat_obj <- load_xenium_seurat(
  rds_path = "path/to/sample.rds",
  celltype_col = "celltype",
  cancer_labels = c("Cancer LumA SC", "Cancer LumB SC")
)

# Automatically detect tumor regions with DBSCAN
seurat_obj <- add_dbscan_polygons(seurat_obj)
auto_polygons <- seurat_obj@images[["fov"]]@boundaries$cancer_regions

# Calculate boundary distances
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, auto_polygons)

# Calculate distance-based separation metrics
results <- calculate_distance_profile(
  seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all"
)

# Visualize
plot_all_metrics(results$profile)
```

## Citation

If you use stGradient in your research, please cite:

> Simons, N.D. (2026). stGradient: R toolkit for distance-based analyses
> of spatial transcriptomics data.

## Contact

For questions, issues, or contributions: - GitHub Issues:
<https://github.com/ndsimons/stGradient/issues> - Email:
<noah.simons@providence.org>
