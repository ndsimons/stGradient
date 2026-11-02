# Plot Distance-Associated Gene Results

Plot Distance-Associated Gene Results

## Usage

``` r
plot_distance_genes(result_list, top_n = 12, sample_label = "")
```

## Arguments

- result_list:

  Output list from \`find_distance_associated_genes\`.

- top_n:

  Maximum number of significant genes shown in trend facets.

- sample_label:

  Optional label for plot titles.

## Value

List with \`volcano\` and \`trends\` ggplot objects. The trends x-axis
extends to the \`max_distance\` used in
\`find_distance_associated_genes()\`.
