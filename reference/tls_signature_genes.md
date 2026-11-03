# Default TLS Signature Genes

Returns the default gene set used to score tertiary lymphoid structures
(TLS) with \[add_tls_score()\]. Genes absent from a panel are dropped at
scoring time.

## Usage

``` r
tls_signature_genes()
```

## Value

Character vector of gene symbols.

## Examples

``` r
tls_signature_genes()
#>  [1] "MS4A1"  "CD19"   "CD22"   "CD3D"   "CD3E"   "CD4"    "CD8A"   "CD8B"  
#>  [9] "CD74"   "CD79A"  "IL7R"   "ITGAE"  "CD1D"   "CD52"   "CD79B"  "FCER2" 
#> [17] "CR2"    "PDCD1"  "PTGDS"  "TRBC2"  "CXCL13" "CXCR5" 
```
