# Package index

## Loading & Setup

Load data and launch the interactive selector

- [`load_xenium_seurat()`](load_xenium_seurat.md) : Load Xenium Seurat
  Object with Stable Cell IDs
- [`launch_stGradient()`](launch_stGradient.md) : Launch stGradient
  Interactive Polygon Selector

## Region Selection

Define tumor regions via DBSCAN or interactive polygons

- [`add_dbscan_polygons()`](add_dbscan_polygons.md) : Automatically
  detect cancer cell regions and create polygons
- [`get_cancer_regions()`](get_cancer_regions.md) : Get Cancer Region
  Polygons
- [`optimize_dbscan_params()`](optimize_dbscan_params.md) : Optimize
  DBSCAN Parameters via Silhouette Score
- [`add_polygon_membership()`](add_polygon_membership.md) : Add Polygon
  Membership to Seurat Object
- [`plot_dbscan_polygons()`](plot_dbscan_polygons.md) : Plot
  DBSCAN-based cancer region polygons

## Distance Metrics

Boundary distances and separation metrics

- [`calculate_boundary_distances()`](calculate_boundary_distances.md) :
  Calculate Boundary Distances
- [`calculate_distance_profile()`](calculate_distance_profile.md) :
  Calculate Distance Profile
- [`calculate_transcriptomic_distance()`](calculate_transcriptomic_distance.md)
  : Calculate Transcriptomic Distance
- [`calculate_composition_difference()`](calculate_composition_difference.md)
  : Calculate Cell Type Composition Difference
- [`calculate_de_strength()`](calculate_de_strength.md) : Calculate
  Differential Expression Strength
- [`find_critical_distance()`](find_critical_distance.md) : Find
  Critical Distance
- [`compute_distance_from_boundary()`](compute_distance_from_boundary.md)
  : Compute Signed Distance From Tumor Boundary

## PCF Analysis

Pair correlation functions outside tumor polygons

- [`build_outside_ppp()`](build_outside_ppp.md) : Build an
  Outside-Polygon Point Pattern for PCF
- [`calculate_pcf_outside_summary()`](calculate_pcf_outside_summary.md)
  : Calculate Outside-Polygon Pair Correlation Summary
- [`plot_pcf_outside_heatmap()`](plot_pcf_outside_heatmap.md) : Plot
  Outside-Polygon PCF Heatmap
- [`plot_pcf_outside_pair_envelopes()`](plot_pcf_outside_pair_envelopes.md)
  : Plot Outside-Polygon PCF Envelopes for All Cell-Type Pairs
- [`plot_pcf_pair_gg()`](plot_pcf_pair_gg.md) : Plot a Single
  Outside-Polygon PCF Pair as a ggplot

## Gene-Distance Analysis

Model gene expression as a function of spatial distance

- [`find_distance_associated_genes()`](find_distance_associated_genes.md)
  : Find Distance-Associated Genes
- [`find_hotspot_distance_genes()`](find_hotspot_distance_genes.md) :
  Find Genes Associated With Distance From T Cell Hotspots
- [`define_tcell_hotspots()`](define_tcell_hotspots.md) : Define
  Intratumoral T Cell Hotspots
- [`plot_distance_genes()`](plot_distance_genes.md) : Plot
  Distance-Associated Gene Results
- [`plot_tcell_hotspots()`](plot_tcell_hotspots.md) : Plot T Cell
  Hotspots on Spatial Coordinates

## Visualization

QC plots, ROI composition, and distance profile visualization

- [`qc_plot_xenium()`](qc_plot_xenium.md) : QC Visualization for Xenium
  Spatial Data
- [`plot_roi_celltype_dist()`](plot_roi_celltype_dist.md) : Plot Cell
  Type Distribution Inside ROI
- [`plot_roi_inout_celltype_stack()`](plot_roi_inout_celltype_stack.md)
  : Plot Stacked Cell Type Composition Inside vs Outside ROI
- [`plot_roi_inout_celltype_dodge()`](plot_roi_inout_celltype_dodge.md)
  : Plot Dodged Cell Type Comparison Inside vs Outside ROI
- [`plot_distance_comparison()`](plot_distance_comparison.md) : Plot
  Distance Profile Comparison
- [`plot_all_metrics()`](plot_all_metrics.md) : Plot All Metrics Faceted
- [`plot_celltype_heatmap()`](plot_celltype_heatmap.md) : Plot Cell Type
  Contributions Heatmap
- [`plot_celltype_proportions()`](plot_celltype_proportions.md) : Plot
  Cell Type Proportions

## Utilities

Internal helpers

- [`parse_svg_path()`](parse_svg_path.md) : Parse SVG Path to
  Coordinates
