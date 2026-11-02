# Plot Cell Type Proportions

Shows proportion trajectories for top cell types

## Usage

``` r
plot_celltype_proportions(
  celltype_contributions,
  top_n = 5,
  span = 0.3,
  level = 0.95,
  ncol = NULL
)
```

## Arguments

- celltype_contributions:

  Data frame from calculate_distance_profile

- top_n:

  Number of top cell types to show

- span:

  Loess smoothing span

- level:

  Confidence level for the loess ribbon

- ncol:

  Number of facet columns. If `NULL` (default), chosen automatically by
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

## Value

ggplot object
