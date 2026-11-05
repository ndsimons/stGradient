# Plot Per-Band PCF Heatmaps

Draws one cell-type by cell-type heatmap of g(r) per distance band from
the output of \[calculate_pcf_by_band()\].

## Usage

``` r
plot_pcf_band_heatmap(
  band_summary,
  distance = NULL,
  value = c("g_value", "g_over_null"),
  midpoint = 1,
  alpha = 0.05
)
```

## Arguments

- band_summary:

  Output of \[calculate_pcf_by_band()\].

- distance:

  Which \`distance\` to plot when several were computed. \`NULL\` uses
  the smallest.

- value:

  \`"g_value"\` (observed g) or \`"g_over_null"\` (observed g divided by
  the null mean). \`"g_over_null"\` removes effects already present in
  the null, e.g. the tumor gradient under \`"label_shuffle"\` with
  \`"anchored"\`. Requires a null model.

- midpoint:

  Midpoint of the diverging color scale.

- alpha:

  Significance level for marking cells with \`\*\` (uses \`p_adj\`). Set
  \`NULL\` to hide markers.

## Value

ggplot object, faceted by band.
