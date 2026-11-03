# Summarize TLS Signature Scores

Sanity check that called TLS score higher than the rest of the tissue.

## Usage

``` r
summarize_tls_score(seurat_obj, score_col = "TLS_score")
```

## Arguments

- seurat_obj:

  Seurat object with \[call_tls()\] and \[add_tls_score()\] results.

- score_col:

  Metadata column with the signature score.

## Value

List with \`by_status\` (in vs. outside TLS), \`by_tls\` (per TLS,
ordered by mean score) and \`wilcox_p\` (two-sided Wilcoxon test
p-value).
