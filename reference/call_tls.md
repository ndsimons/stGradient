# Call Tertiary Lymphoid Structures (TLS)

Clusters lymphoid (B and T) cells with DBSCAN to find candidate immune
aggregates, assigns every cell within \`eps\` of a clustered lymphoid
cell to that aggregate, and calls an aggregate a TLS when it contains
enough co-localized B and T cells. Aggregates containing mature
dendritic cells (e.g. LAMP3+/DC-LAMP+) are flagged as mature-like.

## Usage

``` r
call_tls(
  seurat_obj,
  fov_name = "fov",
  celltype_col,
  b_types = c("B_Cells"),
  t_types = c("CD4+_T_Cells", "CD8+_T_Cells"),
  dc_types = c("LAMP3+_DCs", "IRF7+_DCs"),
  eps = 50,
  minPts = 20,
  min_B = 20,
  min_T = 20,
  min_frac_B = 0.05
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- fov_name:

  Field of view/image name (default: "fov").

- celltype_col:

  Metadata column with cell type labels.

- b_types:

  Cell type labels counted as B cells.

- t_types:

  Cell type labels counted as T cells.

- dc_types:

  Cell type labels counted as mature dendritic cells.

- eps:

  DBSCAN neighborhood radius in coordinate units (microns for Xenium).
  Also the maximum distance for assigning non-lymphoid cells to an
  aggregate.

- minPts:

  DBSCAN minimum lymphoid cells to seed a cluster.

- min_B:

  Minimum B cells for an aggregate to be called a TLS.

- min_T:

  Minimum T cells for an aggregate to be called a TLS.

- min_frac_B:

  Minimum fraction of B cells in the aggregate (excludes pure T cell
  infiltrates).

## Value

Seurat object with \`TLS_id\` (aggregate ID for cells in called TLS,
otherwise \`NA\`) and \`in_TLS\` metadata. The per-aggregate summary is
stored in \`seurat_obj@misc\$tls\[\[fov_name\]\]\`; retrieve it with
\[get_tls_summary()\].

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- call_tls(seurat_obj, celltype_col = "first_type")
get_tls_summary(seurat_obj)
} # }
```
