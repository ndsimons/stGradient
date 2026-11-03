# Load Xenium Seurat Object with Stable Cell IDs

Reads a Seurat object from RDS, enforces unique cell IDs, optionally
computes cancer and ROI metadata columns, and stores coordinates keyed
by cell ID.

## Usage

``` r
load_xenium_seurat(
  rds_path,
  image = NULL,
  celltype_col,
  cancer_labels = NULL,
  roi_cells = NULL,
  roi_col = "inside_polygon",
  store_coords = TRUE,
  coords_name = "tissue_coords"
)
```

## Arguments

- rds_path:

  Path to serialized Seurat object (.rds)

- image:

  Optional image/FOV name to use for coordinates

- celltype_col:

  Metadata column containing cell type annotations

- cancer_labels:

  Optional character vector of labels to mark as cancer

- roi_cells:

  Optional character vector of cell IDs to mark as inside ROI

- roi_col:

  Metadata column name for ROI membership

- store_coords:

  Logical; store processed coordinate data in misc

- coords_name:

  Name used for stored coordinate data in misc

## Value

Seurat object with standardized metadata and optional stored coordinates
