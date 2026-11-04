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
| [TLS Analysis](articles/tls-analysis.html) | Call tertiary lymphoid structures and score a TLS gene signature |
| [Visualization](articles/visualization.html) | QC plots, ROI composition, distance profile plots |

## Quick Start

```r
library(stGradient)

# --- Loading & Setup ---
seurat_obj <- load_xenium_seurat(
  rds_path = "path/to/sample.rds",
  celltype_col = "first_type",
  cancer_labels = c("DCIS 1", "DCIS 2", "Invasive_Tumor", "Prolif_Invasive_Tumor")
)
# Or draw ROIs interactively: launch_stGradient()

# --- Region Selection ---
# Pick DBSCAN parameters by coverage/purity, then build ROI polygons
dbscan_res <- optimize_dbscan_params(seurat_obj, noise_pct = 0.15, b = 0.5, plot = TRUE)
seurat_obj <- add_dbscan_polygons(seurat_obj, eps = dbscan_res$eps[1], minPts = dbscan_res$minPts[1])
auto_polygons <- get_cancer_regions(seurat_obj)
plot_dbscan_polygons(seurat_obj, eps = dbscan_res$eps[1], minPts = dbscan_res$minPts[1])

# --- Distance Metrics ---
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, auto_polygons)
results <- calculate_distance_profile(
  seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "first_type"
)
find_critical_distance(results$profile, "transcriptomic_distance")

# --- PCF Analysis ---
pcf_summary <- calculate_pcf_outside_summary(
  seurat_obj, auto_polygons,
  distance = 200,
  celltype_col = "first_type"
)
plot_pcf_outside_heatmap(pcf_summary, midpoint = 1)

# --- Gene-Distance Analysis ---
seurat_obj <- compute_distance_from_boundary(seurat_obj)
dist_res <- find_distance_associated_genes(seurat_obj, assay = "Xenium", focus = "peritumoral")
dist_plots <- plot_distance_genes(dist_res, top_n = 12)

# --- TLS Analysis ---
seurat_obj <- call_tls(
  seurat_obj,
  celltype_col = "first_type",
  b_types = "B_Cells",
  t_types = c("CD4+_T_Cells", "CD8+_T_Cells")
)
seurat_obj <- add_tls_score(seurat_obj, assay = "Xenium")
plot_tls_score(seurat_obj, highlight = "outline")

# --- Visualization ---
qc_out <- qc_plot_xenium(seurat_obj, roi_col = "inside_polygon", celltype_col = "first_type")
plot_all_metrics(results$profile)
plot_celltype_proportions(results$celltype_contributions, top_n = 5)
```

## Citation

If you use stGradient in your research, please cite:

> Simons, N.D. (2027). stGradient: R toolkit for distance-based analyses of spatial transcriptomics data.

## Contact

For questions, issues, or contributions:
- GitHub Issues: [https://github.com/ndsimons/stGradient/issues](https://github.com/ndsimons/stGradient/issues)
- Email: noah.simons@providence.org
