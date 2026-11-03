# Loading & Setup

## Introduction

This vignette covers the initial steps for working with stGradient:
loading spatial transcriptomics data and launching the interactive
polygon selector.

## Loading Xenium Data

The [`load_xenium_seurat()`](../reference/load_xenium_seurat.md)
function provides a standardized way to load 10x Xenium Seurat objects
with stable cell IDs and cancer/non-cancer annotations.

``` r

library(stGradient)

seurat_obj <- load_xenium_seurat(
  rds_path = "path/to/sample_seuratObject.rds",
  celltype_col = "first_type",
  cancer_labels = c("DCIS 1", "DCIS 2", "Invasive_Tumor", "Prolif_Invasive_Tumor")
)
```

This function:

1.  Reads a Seurat `.rds` file
2.  Creates a stable cell-to-barcode mapping
3.  Adds an `is_cancer` metadata column based on the cell type labels
    you supply

## Next Steps

Once you have your Seurat object loaded, proceed to the [Region
Selection](region-selection.md) vignette to select ROI polygons and
membership assignment.
