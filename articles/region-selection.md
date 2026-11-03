# Region Selection

## Introduction

stGradient provides two approaches for defining regions of interest
(ROIs):

1.  **Interactive selection** via the Shiny app
2.  **Automated detection** using DBSCAN clustering

## Launching the Interactive Polygon Selector

The Shiny app allows you to interactively draw polygons on spatial data
to define ROIs:

``` r

launch_stGradient()
```

### Workflow in the App

1.  **Enter your Seurat object name** — the name of the object in your R
    session
2.  **Enter FOV name** — the field of view name (typically `"fov"`)
3.  **Click “Load Data”** — loads cell coordinates and plots them
4.  **Draw polygons** — use the freehand drawing toolbar on the plot
5.  **Enter polygon object name** — name for the extracted polygon
    object
6.  **Click “Extract Polygons”** — saves polygons to your global
    environment

The resulting polygon object is an `sf` object that can be passed
directly to downstream stGradient functions. ![](stGradient_launch.png)

## DBSCAN-Based Polygon Detection

The [`add_dbscan_polygons()`](../reference/add_dbscan_polygons.md)
function automatically identifies ROIs by clustering cancer-labeled
cells using the DBSCAN algorithm, then wrapping detected clusters in
concave hull polygons.

``` r

library(stGradient)

# One-step: detect polygons and add membership to Seurat object
seurat_obj <- add_dbscan_polygons(seurat_obj)

# Access the polygon object
auto_polygons <- get_cancer_regions(seurat_obj)
```

### Customizing DBSCAN Parameters

You can tune the clustering behavior:

``` r

seurat_obj <- add_dbscan_polygons(
  seurat_obj,
  eps = 35,         # neighborhood radius (microns)
  minPts = 20,      # minimum cluster size
  concavity = 2     # concave hull tightness (lower = tighter)
)
```

### Visualizing Detected Polygons

Use [`plot_dbscan_polygons()`](../reference/plot_dbscan_polygons.md) to
inspect the results before proceeding:

``` r

# Plot with polygons overlaid
plot_dbscan_polygons(seurat_obj)

# Get polygons without plotting (useful for pipelines)
auto_polygons <- plot_dbscan_polygons(seurat_obj, plot = FALSE)
```

![](plot_dbscan_polygons.png)

### Optimizing DBSCAN Parameters

If you’re unsure which `eps` and `minPts` values work best for your
data, use
[`optimize_dbscan_params()`](../reference/optimize_dbscan_params.md) to
perform a grid search evaluated by silhouette score:

``` r

# Search over a grid of eps and minPts values
results <- optimize_dbscan_params(seurat_obj)

# View top parameter combinations
head(results)

# Use the best parameters
best <- results[1, ]
seurat_obj <- add_dbscan_polygons(seurat_obj, eps = best$eps, minPts = best$minPts)
```

You can customize the search ranges:

``` r

results <- optimize_dbscan_params(
  seurat_obj,
  eps_range = seq(20, 80, by = 10),
  minPts_range = c(15, 20, 30, 50),
  max_sample = 30000  # reduce for faster computation
)
```

#### Evaluating Polygon Coverage and Purity

Silhouette score measures cluster separation, but not how well the
resulting polygons capture cancer cells. For candidates with low noise,
build the polygons and score each by `cancer_covered` (fraction of
cancer cells inside a polygon) and `purity` (fraction of cells inside
polygons that are cancer):

``` r

library(ggplot2)
library(tidyr)

# Search over a grid of eps and minPts values
results <- optimize_dbscan_params(seurat_obj)

eval_polys <- function(obj, eps, minPts, fov = "fov") {
  polys <- plot_dbscan_polygons(obj, fov = fov, eps = eps, minPts = minPts, plot = FALSE)
  if (is.null(polys)) return(NULL)
  coords <- Seurat::GetTissueCoordinates(obj, image = fov)
  coords$is_cancer <- as.logical(obj$is_cancer)
  pts    <- sf::st_as_sf(coords, coords = c("x", "y"))
  region <- sf::st_union(sf::st_make_valid(polys))
  inside <- lengths(sf::st_intersects(pts, region)) > 0
  data.frame(eps = eps, minPts = minPts, n_regions = nrow(polys),
             cancer_covered = mean(inside[coords$is_cancer]),
             purity         = mean(coords$is_cancer[inside]))
}

cand <- subset(results, noise_pct < 0.15)[, c("eps", "minPts")]
scores <- do.call(rbind, Map(eval_polys, list(seurat_obj), cand$eps, cand$minPts))
scores[order(-scores$cancer_covered, -scores$purity), ]

scores$label <- paste0(scores$eps, "/", scores$minPts)

b <- 0.5
scores$fbeta <- with(scores, (1 + b^2) * purity * cancer_covered /
                       (b^2 * purity + cancer_covered))

# visualize coverage vs purity
long <- pivot_longer(scores, c(cancer_covered, purity),
                     names_to = "metric", values_to = "value")

ggplot(long, aes(factor(eps), factor(minPts), fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.8) +
  facet_wrap(~ metric, ncol = 1) +
  scale_fill_viridis_c(labels = scales::percent) +
  labs(x = "eps", y = "minPts", fill = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank())
```

`fbeta` combines both metrics into one score; with `b = 0.5`, purity
counts more than coverage.

## Adding Polygon Membership

Once you have polygons (from either the Shiny app or DBSCAN), annotate
each cell as inside or outside:

``` r

# Add inside/outside metadata column
seurat_obj <- add_polygon_membership(seurat_obj, auto_polygons)
```

This adds an `inside_polygon` column to `seurat_obj@meta.data` with
logical values.

## QC Visualization

After assigning membership, use
[`qc_plot_xenium()`](../reference/qc_plot_xenium.md) to visualize the
spatial layout:

``` r

qc_out <- qc_plot_xenium(
  obj = seurat_obj,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  show = "both",
  zoom = TRUE,
  zoom_pad = 0.05,
  grey_outside_roi = TRUE
)

qc_out$p
```

![](qc_out.png)

### ROI Composition Plots

Examine cell type distributions inside vs. outside the region:

``` r

# Distribution inside ROI only
plot_roi_celltype_dist(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "proportion"
)
```

![](plot_roi_celltype_dist.png)

``` r

# Stacked inside vs outside
plot_roi_inout_celltype_stack(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas"
)
```

![](plot_roi_celltype_stack.png)

``` r

#Dodged inside vs outside
plot_roi_inout_celltype_dodge(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "count",
  top_n = 10
)
```

![](plot_roi_inout_celltype_dodge.png)

## Next Steps

With regions defined and membership assigned, proceed to [Distance
Metrics](distance-metrics.md) to calculate boundary distances and
separation metrics.
