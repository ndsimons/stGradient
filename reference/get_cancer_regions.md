# Get Cancer Region Polygons

Retrieves the polygons stored by \[add_dbscan_polygons()\] for a field
of view. Also reads polygons from the FOV \`boundaries\` slot, where
objects created by earlier versions of stGradient stored them.

## Usage

``` r
get_cancer_regions(seurat_obj, fov_name = "fov")
```

## Arguments

- seurat_obj:

  Seurat object processed with \[add_dbscan_polygons()\]

- fov_name:

  Name of the field of view (default: "fov")

## Value

sf object of cancer region polygons, or \`NULL\` with a warning if none
are stored for \`fov_name\`

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- add_dbscan_polygons(seurat_obj)
auto_polygons <- get_cancer_regions(seurat_obj)
} # }
```
