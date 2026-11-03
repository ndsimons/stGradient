# Plot TLS Signature Score Over Space

Colors all cells by TLS signature score (clipped to the 1st-99th
percentile) on a dark tissue background, optionally highlighting called
TLS.

## Usage

``` r
plot_tls_score(
  seurat_obj,
  fov_name = "fov",
  score_col = "TLS_score",
  highlight = c("outline", "ring", "none"),
  label_tls = FALSE,
  concavity = 2,
  point_size = 0.25
)
```

## Arguments

- seurat_obj:

  Seurat object with \[add_tls_score()\] results (and \[call_tls()\]
  results unless \`highlight = "none"\`).

- fov_name:

  Field of view/image name (default: "fov").

- score_col:

  Metadata column with the signature score.

- highlight:

  How to mark called TLS: \`"outline"\` draws a white concave hull
  perimeter around each TLS, \`"ring"\` circles each TLS cell in cyan,
  \`"none"\` shows the score only.

- label_tls:

  Label each TLS with its ID at its centroid.

- concavity:

  Concave hull concavity for \`"outline"\` (lower = tighter).

- point_size:

  Point size for cells.

## Value

ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_tls_score(seurat_obj, highlight = "outline")
} # }
```
