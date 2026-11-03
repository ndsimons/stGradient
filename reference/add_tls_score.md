# Add a TLS Signature Score

Scores each cell for a TLS gene signature with
\[Seurat::AddModuleScore()\], which corrects for background expression
using control gene sets. Genes not in the assay are dropped with a
message.

## Usage

``` r
add_tls_score(
  seurat_obj,
  genes = tls_signature_genes(),
  assay = "Xenium",
  score_name = "TLS_score",
  ctrl = NULL,
  nbin = 24,
  min_genes = 3,
  seed = 1
)
```

## Arguments

- seurat_obj:

  Seurat object.

- genes:

  Signature genes (default: \[tls_signature_genes()\]).

- assay:

  Assay name. If its \`data\` layer is empty,
  \[Seurat::NormalizeData()\] is run with default settings.

- score_name:

  Metadata column to store the score in.

- ctrl:

  Number of control genes per signature gene. Defaults to the number of
  genes present, bounded to 5-50 and capped by the genes available per
  expression bin, which suits small panels.

- nbin:

  Number of expression bins for control gene selection; lower this if
  AddModuleScore fails on small panels.

- min_genes:

  Minimum signature genes that must be present.

- seed:

  Random seed passed to AddModuleScore.

## Value

Seurat object with the score in metadata column \`score_name\`.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- add_tls_score(seurat_obj, assay = "Xenium")
} # }
```
