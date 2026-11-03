# Compare TLS Signature Score In vs. Outside TLS

Plots the per-cell TLS signature score for cells inside called TLS
versus the rest of the tissue, annotated with a two-sided Wilcoxon
rank-sum test p-value.

## Usage

``` r
plot_tls_score_comparison(
  seurat_obj,
  score_col = "TLS_score",
  type = c("beeswarm", "violin")
)
```

## Arguments

- seurat_obj:

  Seurat object with \[call_tls()\] and \[add_tls_score()\] results.

- score_col:

  Metadata column with the signature score.

- type:

  \`"beeswarm"\` (quasirandom points with boxplot; uses ggbeeswarm if
  installed, otherwise jittered points) or \`"violin"\` (violin with
  quartile lines).

## Value

ggplot object.
