# Plot Cell Type Proportions

Shows proportion trajectories for the top cell types, or for cell types
chosen by the user.

## Usage

``` r
plot_celltype_proportions(
  celltype_contributions,
  top_n = 5,
  span = 0.3,
  level = 0.95,
  ncol = NULL,
  celltypes = NULL
)
```

## Arguments

- celltype_contributions:

  Data frame from calculate_distance_profile

- top_n:

  Number of top cell types to show. Ignored when \`celltypes\` is
  supplied.

- span:

  Loess smoothing span

- level:

  Confidence level for the loess ribbon

- ncol:

  Number of facet columns. If \`NULL\` (default), chosen automatically
  by \`ggplot2::facet_wrap()\`.

- celltypes:

  Optional character vector of cell types to plot, shown in the order
  given. If \`NULL\` (default), the \`top_n\` cell types with the
  largest contribution at the peak distance are shown.

## Value

ggplot object
