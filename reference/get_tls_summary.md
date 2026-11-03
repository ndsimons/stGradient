# Get TLS Aggregate Summary

Retrieves the per-aggregate composition table stored by \[call_tls()\].

## Usage

``` r
get_tls_summary(seurat_obj, fov_name = "fov", called_only = TRUE)
```

## Arguments

- seurat_obj:

  Seurat object processed with \[call_tls()\].

- fov_name:

  Field of view/image name (default: "fov").

- called_only:

  If \`TRUE\` (default), return only aggregates called as TLS.

## Value

Data frame with one row per aggregate: cell counts (\`n_total\`,
\`n_B\`, \`n_T\`, \`n_DC\`), fractions, centroid (\`cx\`, \`cy\`),
\`is_TLS\` and \`maturity\`. \`NULL\` with a warning if \[call_tls()\]
has not been run.
