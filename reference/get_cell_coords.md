# Get Cell Coordinates Keyed by Cell ID

Internal helper that returns tissue coordinates matched to Seurat cell
IDs. Uses the Seurat v5 \`cell\` column when present, otherwise the
coordinate rownames. Never matches coordinates to cells by position.

## Usage

``` r
get_cell_coords(seurat_obj, fov_name = "fov")
```

## Arguments

- seurat_obj:

  Seurat object.

- fov_name:

  Field of view/image name.

## Value

Data frame with \`cell_id\`, \`x\`, \`y\` and rownames set to cell IDs,
restricted to cells present in both the FOV and \`seurat_obj\`, in
\`colnames(seurat_obj)\` order.
