# stGradient

<!-- badges: start -->
<!-- badges: end -->
<img src="man/figures/stGradientLogo.png" align="right" width="152" height="175" />

### R toolkit for distance-based analyses of spatial transcriptomics data

stGradient provides tools to interactively select regions of interest (ROIs) and analyze genomic separation at varying distances from tumor boundaries. Designed for spatial transcriptomics analysis, particularly for understanding tumor microenvironment architecture.

## Installation

```r
devtools::install_github("ndsimons/stGradient")
```

## Getting Started

A  [**vignette**](articles/stGradient-workflow.html) is available that covers the core workflow from loading data through visualization.

Additional documentation organized by topic:

| Vignette | Description |
|----------|-------------|
| [Loading & Setup](articles/loading-and-setup.html) | Load Xenium data, launch the interactive Shiny selector |
| [Region Selection](articles/region-selection.html) | Interactive polygon drawing and DBSCAN-based auto-detection |
| [Distance Metrics](articles/distance-metrics.html) | Boundary distances, transcriptomic distance, composition, DE strength |
| [PCF Analysis](articles/pcf-analysis.html) | Pair correlation functions outside tumor polygons |
| [Gene-Distance Analysis](articles/gene-distance-analysis.html) | Model gene expression as a function of distance from boundaries |
| [Visualization](articles/visualization.html) | QC plots, ROI composition, distance profile plots |

## Quick Start

```r
library(stGradient)

# Load a Xenium Seurat object
seurat_obj <- load_xenium_seurat(
  rds_path = "path/to/sample.rds",
  celltype_col = "first_type",
  cancer_labels = c("DCIS 1", "DCIS 2", "Invasive_Tumor", "Prolif_Invasive_Tumor")
)

# Automatically detect ROIs with DBSCAN
seurat_obj <- add_dbscan_polygons(seurat_obj)
auto_polygons <- get_cancer_regions(seurat_obj)

# Calculate boundary distances
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, auto_polygons)

# Calculate distance-based separation metrics
results <- calculate_distance_profile(
  seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "first_type"
)

# Visualize
plot_all_metrics(results$profile)
```

## Citation

If you use stGradient in your research, please cite:

> Simons, N.D. (2027). stGradient: R toolkit for distance-based analyses of spatial transcriptomics data.

## Contact

For questions, issues, or contributions:
- GitHub Issues: [https://github.com/ndsimons/stGradient/issues](https://github.com/ndsimons/stGradient/issues)
- Email: noah.simons@providence.org
