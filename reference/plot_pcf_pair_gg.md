# Plot a Single Outside-Polygon PCF Pair as a ggplot

Plots the pair correlation function (with simulation envelopes) for one
cell-type pair from the list returned by
\[plot_pcf_outside_pair_envelopes()\], as a ggplot object suitable for
embedding in R Markdown documents.

## Usage

``` r
plot_pcf_pair_gg(env_list, celltype1, celltype2)
```

## Arguments

- env_list:

  Named list of envelope objects returned by
  \[plot_pcf_outside_pair_envelopes()\].

- celltype1:

  First cell type label (character).

- celltype2:

  Second cell type label (character). If \`"\<celltype1\> vs
  \<celltype2\>"\` is not in \`env_list\`, the reversed pair is used,
  since unordered envelope lists contain only one direction.

## Value

A ggplot object with the observed g(r) curve and simulation envelope
bands.

## Examples

``` r
if (FALSE) { # \dontrun{
env_list <- plot_pcf_outside_pair_envelopes(seurat_obj, auto_polygons)
plot_pcf_pair_gg(env_list, "TypeA", "TypeB")
} # }
```
