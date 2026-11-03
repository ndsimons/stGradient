# Plot TLS Cellular Composition

Stacked bar chart of cell type fractions within each called TLS.

## Usage

``` r
plot_tls_composition(seurat_obj, celltype_col)
```

## Arguments

- seurat_obj:

  Seurat object processed with \[call_tls()\].

- celltype_col:

  Metadata column with cell type labels.

## Value

ggplot object.
