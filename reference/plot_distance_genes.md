# Plot Distance-Associated Gene Results

Plot Distance-Associated Gene Results

## Usage

``` r
plot_distance_genes(
  result_list,
  top_n = 12,
  sample_label = "",
  max_distance = 200
)
```

## Arguments

- result_list:

  Output list from \`find_distance_associated_genes\`.

- top_n:

  Maximum number of significant genes shown in trend facets.

- sample_label:

  Optional label for plot titles.

- max_distance:

  Upper x-axis limit (distance from boundary) for the trends plot. Set
  to \`NULL\` to show the full distance range.

## Value

List with \`volcano\` and \`trends\` ggplot objects.
