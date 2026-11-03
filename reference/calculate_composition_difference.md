# Calculate Cell Type Composition Difference

Compares cell type distributions between two groups

## Usage

``` r
calculate_composition_difference(
  inside_cells,
  outside_cells,
  seurat_obj,
  group_by
)
```

## Arguments

- inside_cells:

  Character vector of cell IDs inside region

- outside_cells:

  Character vector of cell IDs outside region

- seurat_obj:

  Seurat object

- group_by:

  Metadata column with cell type annotations

## Value

List with total difference and per-celltype breakdown
