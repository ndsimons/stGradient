# Build an Outside-Polygon Point Pattern for PCF

Constructs a marked point pattern (\`ppp\`) using only cells outside the
supplied tumor polygons. This helper is used internally by PCF analysis
functions.

## Usage

``` r
build_outside_ppp(
  seurat_obj,
  polygons,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  fov_name = "fov",
  min_cells = 10,
  window_obj = NULL
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

- window_obj:

  Optional \`spatstat.geom::owin\` window. If \`NULL\`, a bounding-box
  window is built from all coordinates.

## Value

A list with \`ppp_obj\`, \`cells_outside\`, and \`window\`.
