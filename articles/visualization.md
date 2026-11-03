# Visualization

## Introduction

stGradient provides a comprehensive set of plotting functions for
visualizing spatial analysis results. This vignette covers all available
visualizations organized by analysis stage.

## QC and Spatial Overview Plots

### Spatial QC Plot

Visualize cell types and ROI status on the spatial coordinates:

``` r

library(stGradient)

qc_out <- qc_plot_xenium(
  obj = seurat_obj,
  roi_col = "inside_polygon",
  celltype_col = "first_type",
  show = "both",
  zoom = TRUE,
  zoom_pad = 0.05,
  grey_outside_roi = TRUE
)

qc_out$p
```

The `show` parameter accepts `"roi"`, `"celltype"`, or `"both"` for
side-by-side panels.

### ROI Composition Plots

Three complementary views of cell type composition relative to the
region of interest:

``` r

# Cell type distribution inside ROI only
plot_roi_celltype_dist(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "first_type",
  metric = "proportion",
  labels = "both"
)

# Stacked bar: inside vs outside composition
plot_roi_inout_celltype_stack(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "first_type",
  metric = "proportion"
)

# Dodged bar: side-by-side comparison
plot_roi_inout_celltype_dodge(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "first_type",
  metric = "count",
  top_n = 10,
  top_n_by = "inside"
)
```

## Distance Profile Plots

### Transcriptomic Distance Comparison

Compare separation metrics across inside modes:

``` r

plot_distance_comparison(combined_profile)
```

### All Metrics Faceted

View all three separation metrics in a faceted layout:

``` r

plot_all_metrics(combined_profile)
```

### Cell Type Heatmap

Which cell types drive composition differences at each distance:

``` r

plot_celltype_heatmap(celltype_contributions)
```

### Cell Type Proportions

Trajectory of top cell type proportions across distances:

``` r

plot_celltype_proportions(celltype_contributions, top_n = 5)
```

## DBSCAN Polygon Plots

Visualize detected tumor polygons overlaid on cell coordinates:

``` r

plot_dbscan_polygons(seurat_obj)
```

## PCF Plots

### PCF Heatmap

Pairwise cell type co-occurrence outside ROIs:

``` r

plot_pcf_outside_heatmap(pcf_summary, midpoint = 1)
```

### PCF Pair Envelopes

Monte Carlo envelope curves for statistical significance:

``` r

plot_pcf_outside_pair_envelopes(
  seurat_obj = seurat_obj,
  polygons = auto_polygons,
  celltype_col = "first_type",
  nsim = 39,
  file = "pcf_envelopes.png"
)
```

### Single Pair PCF

Focused ggplot2 PCF curve for one cell-type pair:

``` r

plot_pcf_pair_gg(env_list, "TypeA", "TypeB")
```

## Gene-Distance Plots

### Volcano and Trend Plots

``` r

dist_plots <- plot_distance_genes(dist_res, top_n = 12, sample_label = "Sample_01")
dist_plots$volcano
dist_plots$trends
```

## Niche Plots

### Stacked Composition

``` r

plot_niche_stacked_composition(intra$niche_prop_df)
```

### Grouped Bars with Statistics

``` r

plot_niche_grouped_bars(stat_intra, intra_X, intra_y, intra_niches, title = "Intratumoral")
```

### Logistic Regression Coefficients

``` r

plot_niche_coefficients(lr_intra, title = "Intratumoral")
```

## Customization

All plotting functions return `ggplot2` objects that can be further
customized:

``` r

library(ggplot2)

p <- plot_all_metrics(combined_profile)
p + theme_bw() + labs(title = "Custom Title")
```

Use `patchwork` or `cowplot` to combine multiple plots:

``` r

library(patchwork)

p1 <- plot_distance_comparison(combined_profile)
p2 <- plot_celltype_heatmap(celltype_contributions)

p1 / p2
```
