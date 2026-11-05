# Build the Outside-Tumor Observation Window for PCF

Builds the region where outside-tumor cells could have been observed:
the tissue outline minus the tumor polygons. Leaving tumor or empty
background in the window underestimates intensity and biases g(r)
upward.

## Usage

``` r
build_outside_window(
  coords_df,
  polygons,
  tissue_window = c("concave", "convex", "bbox"),
  window_obj = NULL
)
```

## Arguments

- coords_df:

  Data frame with \`x\` and \`y\` for all cells (tumor included).

- polygons:

  sf polygon object representing tumor regions.

- tissue_window:

  Tissue outline used when \`window_obj\` is \`NULL\`: \`"concave"\`
  (concave hull of all cells), \`"convex"\` (convex hull), or \`"bbox"\`
  (bounding box).

- window_obj:

  Optional \`spatstat.geom::owin\` tissue window. Tumor polygons are
  subtracted from it.

## Value

A \`spatstat.geom::owin\` window.
