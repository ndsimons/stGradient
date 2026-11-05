# Plot Outside-Polygon PCF Envelopes for All Cell-Type Pairs

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
  nsim,
  bw = 10,
  window_obj = NULL,
  file = NULL,
  width = 7,
  height = 5,
  res = 150,
  verbose = TRUE,
  tissue_window = c("concave", "convex", "bbox"),
  null_model = c("csr", "band_permutation"),
  band_width = 25,
  envelope_type = c("global", "pointwise"),
  r_range = NULL
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

  Number of simulations for envelopes (required). Global envelopes have
  significance level \`1 / (nsim + 1)\` (0.05 for \`nsim = 19\`, 0.01
  for \`nsim = 99\`); pointwise envelopes have \`2 / (nsim + 1)\` at
  each r.

- bw:

  Kernel bandwidth (microns) passed to \`pcfcross\`. Use the same value
  as in \[calculate_pcf_outside_summary()\] so the envelope curves match
  the heatmap values. \`NULL\` uses the \`spatstat\` default.

- window_obj:

  Optional \`spatstat.geom::owin\` tissue window; tumor polygons are
  subtracted from it.

- file:

  Optional base PNG file path (e.g. \`"pcf_envelopes.png"\`). One PNG is
  written per pair, with the pair name appended to the file stem.

- width, height:

  PNG size in inches.

- res:

  PNG resolution in pixels per inch.

- verbose:

  Logical; print progress messages.

- tissue_window:

  Tissue outline used when \`window_obj\` is \`NULL\`: \`"concave"\`
  (default; concave hull of all cells), \`"convex"\`, or \`"bbox"\`.
  Tumor polygons are always removed from the window.

- null_model:

  Null model for the envelopes: \`"csr"\` or \`"band_permutation"\`. See
  Details.

- band_width:

  Band width (microns) for \`null_model = "band_permutation"\`.

- envelope_type:

  \`"global"\` (default) or \`"pointwise"\`. See Details.

- r_range:

  Distance range \`c(lower, upper)\` in microns. \`upper\` sets the
  largest distance computed (default: the \`spatstat\` default maximum).
  For global envelopes the test covers \`\[lower, upper\]\`; \`lower\`
  defaults to \`bw\` (or 5

Named list of envelope objects. Each has attributes \`null_model\`,
\`envelope_type\`, \`alpha\` (significance level), and \`r_range\`;
global envelopes also have \`p_value\` (MAD test). Computes and plots
envelope-based pair correlation curves for all requested cell-type pairs
outside polygons. When \`file\` is provided, each pair is saved as its
own PNG named \`\<file stem\>\_\<type1\>\_vs\_\<type2\>.png\`. Two null
models are available: \* \`"csr"\`: complete spatial randomness. Points
are placed uniformly at random in the outside-tumor window with
independent labels. This tests clustering of the cells themselves as
well as co-localization, and tissue is rarely uniform. \*
\`"band_permutation"\`: cell positions are kept and cell-type labels are
permuted among cells in the same distance band from the tumor boundary
(bands of width \`band_width\`). This tests whether two types
co-localize beyond what their shared distance-to-tumor gradients
explain. Use \`band_width = Inf\` for plain random labelling across all
outside cells.Two envelope types are available: \* \`"global"\`
(default): a simultaneous envelope (\`spatstat\` \`global = TRUE\`,
maximum absolute deviation). The observed curve leaves it anywhere in
\`r_range\` with probability \`1 / (nsim + 1)\` under the null, so it
can be read across all distances at once without multiple testing. The
band has the same width at every distance. A matching MAD test p-value,
computed from the same simulations, is stored in the \`p_value\`
attribute. With \`null_model = "band_permutation"\`, \`2 \* nsim\`
permutations are run: the second set estimates the null mean curve. \*
\`"pointwise"\`: separate envelopes at each distance, with a two-sided
significance level of \`2 / (nsim + 1)\` at each single r. Checking many
distances inflates the false-positive rate.
