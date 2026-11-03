# Plot Outside-Polygon PCF Envelopes for All Cell-Type Pairs

Computes and plots envelope-based pair correlation curves for all
requested cell-type pairs outside polygons. When \`file\` is provided,
each pair is saved as its own PNG named \`\<file
stem\>\_\<type1\>\_vs\_\<type2\>.png\`.

## Usage

``` r
plot_pcf_outside_pair_envelopes(
  seurat_obj,
  polygons,
  celltype_col,
  fov_name = "fov",
  min_cells = 10,
  include_self_pairs = FALSE,
  ordered_pairs = FALSE,
  nsim = 10,
  bw = 10,
  window_obj = NULL,
  file = NULL,
  width = 7,
  height = 5,
  res = 150,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- polygons:

  sf polygon object representing ROIs.

- celltype_col:

  Metadata column containing cell type labels.

- fov_name:

  Field of view/image name passed to \`Seurat::GetTissueCoordinates\`.

- min_cells:

  Minimum number of cells required per cell type.

- include_self_pairs:

  Logical; include same-type pairs.

- ordered_pairs:

  Logical; if \`TRUE\`, computes ordered pairs.

- nsim:

  Number of simulations for envelopes.

- bw:

  Bandwidth passed to \`pcfcross\` in envelope evaluation.

- window_obj:

  Optional \`spatstat.geom::owin\` window.

- file:

  Optional base PNG file path (e.g. \`"pcf_envelopes.png"\`). One PNG is
  written per pair, with the pair name appended to the file stem.

- width, height:

  PNG size in inches.

- res:

  PNG resolution in pixels per inch.

- verbose:

  Logical; print progress messages.

## Value

Named list of envelope objects.
