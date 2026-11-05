# Plot PCF Across Distance Bands

Plots g(r) against distance from the tumor (band midpoint) for each
cell-type pair, with the null 95

## Usage

``` r
plot_pcf_band_profile(
  band_summary,
  pairs = NULL,
  distance = NULL,
  alpha = 0.05
)
```

## Arguments

- band_summary:

  Output of \[calculate_pcf_by_band()\].

- pairs:

  Optional character vector of pairs to show, written \`"TypeA vs
  TypeB"\`. \`NULL\` shows all pairs.

- distance:

  Which \`distance\` to plot when several were computed. \`NULL\` uses
  the smallest.

- alpha:

  Significance level for filled points (uses \`p_adj\`).

## Value

ggplot object, faceted by pair.
