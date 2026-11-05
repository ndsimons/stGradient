<!--
Target: Bioinformatics, Original Paper (suggested section: Gene expression).
Placeholders to fill are written as [[NAME]]; search for "[[" to find them all.
Items marked [confirm] are inferred from the code/figures and should be checked.
-->

# stGradient: an R toolkit for distance-resolved analysis of the tumor–stroma interface in single-cell spatial transcriptomics

Noah D. Simons^1,\*^, Anton Kravchenko^1^, Alexa Dowdell^1^

^1^ [[AFFILIATION]], Providence [[CITY, STATE, COUNTRY]]

\*To whom correspondence should be addressed: noah.simons@providence.org

---

## Abstract

**Motivation:** Imaging-based spatial transcriptomics platforms such as 10x Genomics Xenium measure gene expression in hundreds of thousands of cells while keeping their tissue coordinates. Many questions in tumor immunology concern the interface between malignant epithelium and the surrounding stroma. How wide is the zone the tumor influences? Which cells line the boundary? Do immune cells reach the tumor, or are they held in the stroma? Answering these questions requires a reproducible definition of the tumor boundary and analyses organized around distance from it. Current general-purpose spatial toolkits provide neighborhood and domain analyses but no integrated, boundary-centric workflow.

**Results:** We present stGradient, an R package built on Seurat that (i) defines tumor regions either interactively or automatically, using DBSCAN clustering with concave hulls and a coverage–purity parameter search; (ii) assigns every cell a signed distance to the nearest tumor boundary; and (iii) uses this coordinate to quantify compositional and transcriptional separation as a function of distance, cell–cell spatial association outside the tumor, distance-associated genes, and tertiary lymphoid structures. Applied to a Xenium breast cancer section, stGradient showed a macrophage-enriched peritumoral cuff within about 50 µm of tumor nests. T cells were enriched farther out, and CD8^+^ T cells co-localized strongly with these macrophages. This pattern is consistent with T cells being retained in the stroma by macrophage interactions.

**Availability and implementation:** stGradient is freely available under the MIT license at https://github.com/ndsimons/stGradient, with documentation and vignettes at [[PKGDOWN_URL]].

**Contact:** noah.simons@providence.org

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

## 1 Introduction

A solid tumor is not a homogeneous mass of malignant cells. It is an organ-like structure whose behavior depends heavily on the interface between cancer epithelium and the surrounding stroma. Much of what determines whether an anti-tumor immune response succeeds happens at or near this interface. T cells must leave the vasculature in the stroma, cross a margin rich in fibroblasts, myeloid cells and extracellular matrix, and physically contact tumor cells. Tumors are therefore classified as immune-inflamed, immune-excluded or immune-desert (Chen and Mellman, 2017). The immune-excluded phenotype, in which T cells accumulate in the stroma but fail to infiltrate tumor nests, is associated with resistance to immune checkpoint blockade (Mariathasan *et al.*, 2018). Exclusion is a spatial phenotype by definition. Bulk or dissociated single-cell sequencing cannot detect it, and in multiplexed imaging it is usually scored with coarse compartment labels (Keren *et al.*, 2018).

Several mechanisms for exclusion have been proposed. They include physical barriers formed by cancer-associated fibroblasts and dense matrix, TGF-β signaling (Mariathasan *et al.*, 2018), and long-lasting interactions between CD8^+^ T cells and tumor-associated macrophages in the stroma. These interactions slow T-cell migration and keep T cells from reaching tumor cells (Peranzoni *et al.*, 2018). Distinguishing these mechanisms in patient tissue requires single-cell measurements that are explicitly organized around distance from the tumor boundary.

Imaging-based spatial transcriptomics platforms now provide such data. Platforms such as 10x Genomics Xenium resolve hundreds of transcripts in hundreds of thousands of segmented cells per section (Janesick *et al.*, 2023; Moses and Pachter, 2022). General-purpose frameworks such as Seurat (Hao *et al.*, 2024), Giotto (Dries *et al.*, 2021), Squidpy (Palla *et al.*, 2022) and the Bioconductor SpatialExperiment ecosystem (Righelli *et al.*, 2022) provide infrastructure for these data. They also offer neighborhood enrichment, spatially variable gene detection and tissue domain segmentation (e.g. Singhal *et al.*, 2024). Mature point-process libraries such as spatstat (Baddeley *et al.*, 2015) provide rigorous spatial statistics. Studies of tumor–microenvironment interfaces have shown that boundary-adjacent cells have distinct transcriptional states (Hunter *et al.*, 2021). However, analyses of this kind are usually written ad hoc. Each study defines the boundary, the distance coordinate and the inside/outside contrast differently, which makes results hard to reproduce or compare across samples.

Three practical problems recur. First, the boundary itself must be defined. In the "inside/outside" framing, the boundary is a model imposed on the data, not a property of it. A boundary that is too loose folds stroma into the tumor and dilutes every contrast. A boundary that is too tight mislabels edge cancer cells as stroma and creates spurious peritumoral signal. Second, the spatial scale of the interface is usually unknown in advance, yet analyses typically fix it at an arbitrary cut-off such as 50 or 100 µm. Third, a single interface question has several facets: who is present, how cells are arranged, which genes change and whether organized immune structures are present. Each facet calls for a different statistical tool, and these tools rarely share a coordinate system.

stGradient addresses these problems with a single organizing variable: the signed Euclidean distance from every cell to the nearest tumor boundary. It is negative inside tumor nests and positive outside, in microns. The package provides (i) interactive and automated tumor region definition with a transparent parameter-selection criterion; (ii) distance profiles of compositional and transcriptional separation, together with a data-driven estimate of interface width; (iii) pair correlation function (PCF) analysis restricted to the extratumoral compartment; (iv) per-gene models of expression as a function of boundary distance; and (v) calling and scoring of tertiary lymphoid structures (TLS). All functions operate on standard Seurat objects and sf geometries, so they slot into existing single-cell workflows. Here we describe the methods and their implementation. We then show them on a public Xenium breast cancer section, where the combination of analyses points to a macrophage-rich peritumoral cuff that may retain CD8^+^ T cells in the stroma.

## 2 Methods

### 2.1 Overview and data model

stGradient is implemented in R (≥ [[R_VERSION]]) and takes a Seurat object (v4 or v5) with a spatial field of view (FOV) as input. Cell centroids are retrieved with `Seurat::GetTissueCoordinates()` and indexed by cell identifier. All per-cell results are written back to the object's metadata in `colnames()` order, with `NA` for cells lacking coordinates. Region polygons are represented as `sf` simple-feature collections (Pebesma, 2018). Automatically derived polygons are stored in the object's `@misc` slot, keyed by FOV, so that they are carried along with the data through subsetting and saving. The helper `load_xenium_seurat()` reads a serialized object and enforces unique cell identifiers. It also derives a logical `is_cancer` column from a user-supplied set of malignant cell-type labels and caches coordinates. Every downstream step depends only on this column, a cell-type column and the expression matrix, so stGradient is agnostic to how cells were segmented and annotated. Table 1 lists the exported functions by analysis module, and Figure 1 shows the region-definition workflow.

### 2.2 Tumor region definition

**Manual delineation.** `launch_stGradient()` opens a Shiny application built on plotly (Sievert, 2020). The application plots every cell colored by `is_cancer` and lets the user draw free-form closed paths, rectangles or circles (Fig. 1A). Drawn shapes are parsed from plotly's SVG path representation into closed polygons. Circles are discretized into 64 vertices. The shapes are returned to the R session as an `sf` object. Manual delineation is appropriate when a pathologist's judgment should define the region of interest (ROI), for example a specific invasive front, or when a necrotic area should be excluded. It can also be used to correct automated boundaries.

**Automated delineation.** `add_dbscan_polygons()` clusters the coordinates of cancer cells with DBSCAN (Ester *et al.*, 1996), as implemented in the dbscan package (Hahsler *et al.*, 2019). A cancer cell is a core point if at least `minPts` cancer cells lie within radius `eps`. Clusters are chains of density-connected core points. Isolated cancer cells are labeled noise and discarded. This is intended: single infiltrating or misannotated cells should not define a tumor nest. A concave hull is then fitted to each cluster with the concaveman algorithm (Park and Oh, 2012), producing one polygon per nest (Fig. 1B). The `concavity` parameter controls how tightly the hull follows the cells. Because each nest receives its own polygon, stroma between nests is correctly classed as extratumoral. A single hull around all cancer cells would instead label it intratumoral.

**Parameter selection.** DBSCAN results depend strongly on `eps` and `minPts`. `optimize_dbscan_params()` performs a grid search over both parameters. For every combination it records the number of clusters, the fraction of cancer cells labeled noise and the mean silhouette width of the clustered coordinates (Rousseeuw, 1987). Silhouette width is computed on a random subsample (default 46,000 cells) to bound memory. Silhouette is a purely geometric criterion that rewards compact, well-separated clusters, but it does not reflect whether polygons follow the histological tumor edge. stGradient therefore scores each combination with two region-level quantities (Fig. 1C). *Coverage* is the fraction of all cancer cells that fall inside any polygon. *Purity* is the fraction of all cells inside polygons that are cancer cells. Only combinations whose noise fraction is below a user threshold (default 15%) are scored. Coverage and purity are combined into an F~β~ score,

F~β~ = (1 + β^2^) · purity · coverage / (β^2^ · purity + coverage),

and the grid is ranked by F~β~, with silhouette width as a tie-breaker. The default β = 0.5 weights purity more heavily than coverage. Impure polygons contaminate the "inside" reference that every subsequent analysis relies on, whereas modest under-coverage only removes a few edge cells from the analysis. The search returns the full table and a heatmap of coverage, purity and F~β~ with the selected combination outlined (Fig. 1D).

### 2.3 Signed boundary distance

Given a polygon set *P*, `calculate_boundary_distances()` converts the polygons to their boundary linestrings and computes, for each cell *i*, the Euclidean distance *d~i~* to the nearest polygon edge with `sf::st_distance()`. It then tests whether the cell intersects any polygon. The signed distance is *s~i~* = −*d~i~* for cells inside a polygon and *s~i~* = +*d~i~* otherwise. `compute_distance_from_boundary()` stores the same quantity (`dist_to_boundary`) together with its absolute value, and by default reads the polygons saved by `add_dbscan_polygons()`. This single coordinate underlies every analysis below.

### 2.4 Distance profiles and interface width

The central question of the package is how far the tumor's influence extends into the surrounding tissue. `calculate_distance_profile()` evaluates a grid of thresholds *d* (default 10–200 µm in 10 µm steps). At each threshold it defines an outside group of all cells with 0 < *s* ≤ *d* and an inside group selected in one of two modes. In `"all"` mode the inside group is every intratumoral cell (*s* < 0), which serves as a fixed reference. In `"distance"` mode it is only intratumoral cells within *d* of the boundary (−*d* ≤ *s* < 0), giving a symmetric band of width *d* on each side of the edge. Thresholds with fewer than 10 cells in either group are skipped. Three complementary separation metrics are computed for each pair of groups.

*Transcriptomic distance*, *D*~tx~ = 1 − cor(x̄~in~, x̄~out~), is one minus the Pearson correlation between the pseudobulk mean log-normalized expression profiles of the two groups. It is a global, annotation-free summary of how different the two sides look transcriptionally.

*Composition difference*, *D*~comp~ = Σ~k~ (p~k~^in^ − p~k~^out^)^2^, is the squared Euclidean distance between the cell-type proportion vectors of the two groups. Each term is returned as that cell type's contribution, together with the inside and outside proportions at every *d*. This shows which populations drive separation at each distance.

*Differential expression (DE) strength* is the mean |log~2~FC| of the 100 genes with the largest absolute fold change from `Seurat::FindMarkers()`, run as a Wilcoxon test with min.pct = 0.1 and no fold-change threshold. Unlike *D*~tx~, which can wash out a small number of strongly changing genes, DE strength is sensitive to gene-restricted programs such as chemokines or matrix genes.

The outside group is a cumulative shell, so each threshold adds cells to the previous set and the curves tend to rise monotonically before saturating. The value of a metric at threshold *d* is therefore a weighted average over the shell (0, *d*]. For a cell type whose local density peaks near the boundary and then falls, the cumulative proportion rises and then flattens or declines. For a cell type whose density keeps increasing with distance, the cumulative proportion keeps rising. `find_critical_distance()` summarizes a profile by the threshold at which the metric is maximal and by the *plateau distance*, the smallest *d* at which the metric reaches 90% of its maximum. Because the cumulative design makes the maximum typically fall at the largest threshold tested, the plateau distance is the more informative quantity. It serves as an operational, data-driven estimate of the width of the tumor-influenced zone, and can be used to define "peritumoral" in downstream analyses instead of an arbitrary cut-off.

### 2.5 Spatial association outside the tumor

Composition describes who is present but not how cells are arranged. `calculate_pcf_outside_summary()` builds a multitype point pattern (spatstat `ppp`) from the coordinates of cells lying outside all tumor polygons, marked by cell type. Types with fewer than `min_cells` cells (default 10) are dropped. For every pair of cell types (*i*, *j*), optionally including self-pairs, it estimates the cross-type pair correlation function *g*~ij~(*r*) with `spatstat.explore::pcfcross()` and Ripley isotropic edge correction, then reports *g*~ij~ at a user-chosen distance *r*. The PCF is the density of type-*j* cells at distance *r* from a typical type-*i* cell, divided by the density expected if the two types were placed independently. Values above 1 indicate co-localization and values below 1 indicate segregation. Unlike Ripley's *K*, the PCF is non-cumulative, so it localizes associations to a specific spatial scale. `plot_pcf_outside_heatmap()` displays the resulting matrix with a diverging scale centered on 1. `plot_pcf_outside_pair_envelopes()` plots full *g*~ij~(*r*) curves against pointwise simulation envelopes under complete spatial randomness (CSR), computed with `spatstat.explore::envelope()` with a user-defined number of simulations and smoothing bandwidth (default 10 µm). Restricting the analysis to extratumoral cells keeps dense cancer nests from dominating the pattern. It also targets the question of how immune and stromal cells are organized in the space they actually occupy.

### 2.6 Distance-associated genes

`find_distance_associated_genes()` identifies genes whose expression changes continuously with distance from the boundary. Cells are selected according to `focus`: `"peritumoral"` (*s* > 0), `"intratumoral"` (*s* ≤ 0, using depth |*s*|) or `"both"` (signed distance across the boundary). The distance range is trimmed at the smaller of a quantile (default 95th) and an absolute cap (default 200 µm), so that sparse, far-field cells do not dominate the fit. Genes detected in fewer than `min_cells_expressing` cells (default 10) are removed. If the normalized data layer is empty, `Seurat::NormalizeData()` is run first. For each gene, an ordinary least-squares model of log-normalized expression on distance is fitted. The function records the slope (change per µm), its *t* statistic and *P* value, and *R*^2^. Spearman's ρ is also computed independently as a rank-based, outlier-robust measure of monotonic association. *P* values from both tests are adjusted with the Benjamini–Hochberg procedure. For up to 50 significant genes, expression is averaged within quantile-based distance bins (default 10) to produce trend profiles. `plot_distance_genes()` returns a volcano plot of slope against adjusted significance and faceted trend plots for the top genes.

### 2.7 Tertiary lymphoid structures

TLS are ectopic lymphoid aggregates of B cells, T cells and dendritic cells. Their presence and maturity are associated with survival and with response to immune checkpoint blockade (Sautès-Fridman *et al.*, 2019; Cabrita *et al.*, 2020; Helmink *et al.*, 2020; Petitprez *et al.*, 2020). `call_tls()` clusters B and T cells (user-defined label sets) with DBSCAN (default `eps` = 50 µm, `minPts` = 20) to find candidate lymphoid aggregates. It then assigns every cell, of any type, that lies within `eps` of a clustered lymphoid cell to that aggregate, so that embedded stromal, myeloid and dendritic cells are counted. An aggregate is called a TLS when it contains at least `min_B` B cells and `min_T` T cells (default 20 each) and B cells make up at least `min_frac_B` (default 5%) of its cells. The B-cell requirement distinguishes TLS from diffuse T-cell infiltrates. TLS containing mature dendritic cells (by default LAMP3^+^ or IRF7^+^ DCs) are flagged as "mature-like". Per-aggregate cell counts, fractions, centroids and approximate areas are stored with the object, and cells are labeled with `TLS_id` and `in_TLS`.

Independently of the structural call, `add_tls_score()` scores every cell for a 22-gene TLS signature with `Seurat::AddModuleScore()`, which subtracts the mean expression of expression-matched control genes. The signature covers B-cell, T-cell, follicular/germinal-center and antigen-presentation genes (`tls_signature_genes()`). Genes absent from a targeted panel are dropped with a message. The number of control genes is bounded automatically to suit small panels. `summarize_tls_score()` and `plot_tls_score_comparison()` compare scores inside and outside called TLS. `plot_tls_map()`, `plot_tls_composition()` and `plot_tls_score()` visualize calls, their composition and the score on the tissue.

### 2.8 Visualization and implementation

All plotting functions return ggplot2 objects (Wickham, 2016) that can be modified further. They cover spatial QC of ROI membership and cell types (`qc_plot_xenium()`), inside/outside composition (`plot_roi_celltype_dist()`, `plot_roi_inout_celltype_stack()`, `plot_roi_inout_celltype_dodge()`), distance profiles (`plot_all_metrics()`, `plot_distance_comparison()`) and per-cell-type contributions and proportions (`plot_celltype_heatmap()`, `plot_celltype_proportions()`). Heavy optional dependencies (spatstat, cluster, FNN) are listed under Suggests and checked at run time. The package includes eight vignettes, which cover loading, region selection, distance metrics, PCF analysis, gene–distance analysis, TLS analysis and visualization, together with an end-to-end workflow. A parameterized R Markdown file regenerates every figure from a user-supplied object. [[TEST_SUITE_SENTENCE — e.g. "Unit tests (testthat) cover ...", once added]]

### 2.9 Data and analysis settings

We analyzed the publicly available 10x Genomics Xenium FFPE human breast cancer dataset (Janesick *et al.*, 2023), [[SAMPLE / REPLICATE]] [confirm], comprising [[N_CELLS]] cells and [[N_GENES]] genes [confirm: 313-gene panel]. Cell types were annotated by [[ANNOTATION_METHOD]] [confirm: RCTD (Cable *et al.*, 2022), `first_type` labels] against [[REFERENCE]]. Cells labeled DCIS 1, DCIS 2, Invasive_Tumor or Prolif_Invasive_Tumor were treated as malignant. DBSCAN parameters were searched over `eps` ∈ {15, 20, …, 50} µm and `minPts` ∈ {5, 10, …, 45}, with a 15% noise ceiling and β = 0.5. Tumor polygons for downstream analysis were built with `eps` = [[EPS_USED]] µm, `minPts` = [[MINPTS_USED]] and concavity = 2. Distance profiles used thresholds of 10–200 µm in both inside modes. PCF heatmaps were evaluated at *r* = [[PCF_R]] µm [confirm: 200 µm], and pairwise envelopes used [[NSIM]] CSR simulations [confirm: 39] with bandwidth 10 µm. Gene–distance analysis used the peritumoral focus with a 200 µm cap. TLS calling used the default parameters above. Analyses were run on [[HARDWARE]] with R [[R_VERSION]], Seurat [[SEURAT_VERSION]] and stGradient 0.1.0.

## 3 Results

### 3.1 Transparent, tunable definition of tumor nests

The section contains a mixture of large invasive tumor regions and discrete DCIS ducts embedded in stroma (Fig. 1B, top). The coverage–purity search exposed a clear trade-off (Fig. 1D). Small `eps` with large `minPts` produced high purity (up to 0.98) but low coverage (< 0.5), because many edge and scattered cancer cells were discarded as noise. Large `eps` with small `minPts` raised coverage to about 0.95 but lowered purity to about 0.72–0.75, because hulls bridged neighboring nests and swallowed intervening stroma, the failure mode sketched in Fig. 1C. The F~0.5~ criterion selected `eps` = 25 µm and `minPts` = 10 (coverage 0.83, purity 0.88, F~0.5~ = 0.87). A broad ridge of neighboring settings scored within a few hundredths of this optimum, suggesting that downstream results are not sensitive to the exact choice. [[SENSITIVITY_SENTENCE — e.g. "Plateau distances varied by ≤ X µm across the top five settings (Supplementary Table S1)."]] The resulting [[N_POLYGONS]] polygons outlined individual DCIS ducts and invasive nests separately (Fig. 1B, bottom). Stroma between nests was therefore classified as extratumoral. The inside compartment was dominated by malignant cell types, while stromal cells, T cells, endothelial cells, myoepithelial cells and macrophages were found predominantly outside (Fig. 2A). The interactive selector (Fig. 1A) produced polygons that are interchangeable with the automated ones in all downstream functions.

### 3.2 The tumor–stroma interface is compositionally sharp but transcriptionally graded

Distance profiles separated two kinds of change at the interface (Fig. 2B). Composition difference and DE strength rose steeply over the first 20–30 µm and then saturated. Their plateau distances were [[PLATEAU_COMP]] µm and [[PLATEAU_DE]] µm respectively in `"all"` mode. Transcriptomic distance, in contrast, kept rising across the full 200 µm range (plateau [[PLATEAU_TX]] µm). The tumor therefore marks a sharp change in which cells are present within a few cell diameters of its edge. The pseudobulk expression profile of the surrounding tissue, however, continues to diverge from the tumor over a much broader zone. Comparing inside modes showed that composition and DE strength were consistently higher when the boundary band (`"distance"` mode) rather than the whole tumor was used as the reference. The tumor edge therefore differs from the tumor core in composition and in its most strongly changing genes. Transcriptomic distance was nearly identical in both modes.

Per-cell-type contributions showed what drives the compositional separation (Fig. 2C). Over most of the range, the loss of invasive tumor cells and the gain of stromal cells dominated. ACTA2^+^ myoepithelial cells contributed most at 10–20 µm, as expected for the myoepithelial layer that surrounds DCIS ducts. This serves as an internal check that the polygons follow true ductal boundaries.

The cumulative proportion trajectories (Fig. 2D) showed a distinct spatial ordering of immune populations. Macrophages_1 rose from about 3% inside the tumor to a maximum of about 9.5% at roughly 40–60 µm and then flattened or declined slightly. Because the outside group is cumulative, this shape means that Macrophages_1 density is highest immediately adjacent to the boundary. In contrast, CD4^+^ and CD8^+^ T-cell proportions were near zero inside the tumor (< 1%) and rose steadily to the edge of the analyzed range (≈10% and ≈6% at 200 µm), with no plateau. T-cell density therefore keeps increasing with distance from the tumor, beyond the macrophage-rich zone. Macrophages_2 followed the T-cell pattern rather than the Macrophages_1 pattern. The tissue around tumor nests is thus organized as layers: a thin myoepithelial rim, then a Macrophages_1-enriched peritumoral cuff within about 50 µm, then a T-cell-rich stroma farther out. The T cells are largely absent from the tumor itself.

### 3.3 CD8^+^ T cells co-localize with peritumoral macrophages

The layered arrangement raises the question of whether T cells and macrophages are simply adjacent populations or whether they interact. Outside-tumor PCF analysis (Fig. 3A) showed a coherent immune block: B cells, CD4^+^ and CD8^+^ T cells, and LAMP3^+^ and IRF7^+^ dendritic cells were all mutually co-localized, consistent with organized lymphoid aggregates (Section 3.5). The strongest self-associations involved malignant and myoepithelial types that fall outside the polygons (DCIS 1, DCIS 2, KRT15^+^ myoepithelial cells), which mark small ducts and isolated nests rejected as DBSCAN noise. Macrophages_1 were moderately co-localized with both CD4^+^ and CD8^+^ T cells at the summary scale.

The full CD8^+^ T cell–Macrophages_1 PCF (Fig. 3B) showed a strong, short-range association. *g*(*r*) reached about 10 at the shortest distances, remained around 4 at 25–30 µm (≈ 1–2 cell diameters, the scale of direct cell contact), and decayed slowly toward 1 over several hundred microns. It lay far above the CSR envelope at every distance up to about 1 mm. A key alternative explanation is shared habitat: because both populations are confined to stroma, they would appear co-localized even with no specific interaction. To address this, we used CD8^+^ T cells and generic stromal cells, which occupy the same extratumoral space, as a reference pair (Fig. 3C). That PCF fell to about 1.8 by 20 µm and stayed below 2 thereafter. Over the contact range of 20–50 µm, CD8^+^ T cells were therefore found about twice as often near Macrophages_1 as near stromal cells, relative to their respective densities. [[RANDOM_LABEL_SENTENCE — e.g. "Under a random-labelling null that conditions on cell positions, the CD8–Macrophages_1 association remained significant (pointwise P < 0.05 at r ≤ X µm)."]]

Together, Sections 3.2 and 3.3 support the following hypothesis. Peritumoral Macrophages_1 form a cuff around tumor nests, and CD8^+^ T cells engage these macrophages in the stroma instead of reaching tumor cells. In this way, the T cells are effectively held in the extratumoral compartment. This mirrors the macrophage-mediated T-cell exclusion described by live imaging in human and murine tumors (Peranzoni *et al.*, 2018) and fits the immune-excluded phenotype (Chen and Mellman, 2017). stGradient does not establish mechanism or direction of causality from a single static section. However, it turns a qualitative impression from a tissue image into quantitative, distance-resolved statistics that can be compared across patients.

### 3.4 Distance-associated genes reflect the same layered structure

Peritumoral gene–distance analysis identified [[N_SIG_GENES]] of [[N_TESTED_GENES]] tested genes as significantly associated with distance (BH-adjusted *P* < 0.05). For 88 genes the adjusted *P* value fell below machine precision (Fig. 4A). This reflects the very large number of cells, which are not independent replicates. Genes that decreased with distance into the stroma were dominated by luminal epithelial and tumor markers (KRT7, KRT8, EPCAM, FOXA1, GATA3, FASN, CCND1; Fig. 4B). This most likely reflects residual tumor cells and transcript spill-over from segmentation at the interface rather than stromal regulation. stGradient flags these genes for this reason (Section 4). Genes that increased with distance included T-cell and pan-immune markers (CD3E, TRAC, IL7R, PTPRC, CYTIP, CXCR4) and stromal/fibroblast genes (CXCL12, CCDC80, PTGDS). The steady rise of T-cell transcripts over 200 µm independently reproduces the T-cell gradient seen in the composition profiles (Fig. 2D). The co-increase of CXCL12, a chemokine that retains CXCR4^+^ lymphocytes, with CXCR4 itself suggests a candidate chemotactic axis that could contribute to T-cell positioning in the distal stroma. Because these trends conflate changes in composition with changes in cell state, stGradient supports repeating the analysis within a single cell type (e.g. CD8^+^ T cells or Macrophages_1) to test whether the macrophages in the cuff or the T cells near them show distance-dependent state changes. [[WITHIN_TYPE_RESULT — optional, if run]]

### 3.5 Tertiary lymphoid structures are present in the peritumoral stroma

`call_tls()` identified 19 TLS across the section (Fig. 5A). They varied widely in size and shape, from compact aggregates to elongated structures that wrap around ducts. All contained B cells (≈5–20% of cells) together with a larger CD4^+^ T-cell fraction, variable CD8^+^ T-cell, macrophage and stromal fractions, and dendritic cells (Fig. 5B). [[N_MATURE]] of the 19 contained LAMP3^+^ or IRF7^+^ DCs and were flagged as mature-like. The TLS signature score was markedly elevated inside called TLS (median [[MEDIAN_IN]] vs. [[MEDIAN_OUT]]; Wilcoxon *P* < 2 × 10^−16^; Fig. 5D), and its spatial distribution closely matched the structural calls (Fig. 5C). As discussed below, this comparison is a consistency check rather than independent validation. Several high-scoring regions outside called TLS, particularly rings around DCIS ducts, are candidates for smaller or less organized lymphoid infiltrates below the calling thresholds. The TLS were located at a median signed distance of [[TLS_MEDIAN_DIST]] µm from the nearest tumor boundary [[optional]], within the T-cell-rich stroma rather than the macrophage-rich cuff. This places these organized immune niches on the far side of the macrophage layer.

### 3.6 Performance

On the [[N_CELLS]]-cell section, region definition with a [[GRID_SIZE]]-combination parameter search took [[TIME_OPT]], polygon construction and signed distances [[TIME_DIST]], a 20-threshold distance profile with all three metrics [[TIME_PROFILE]], all-pairs PCF summary for [[N_TYPES]] cell types [[TIME_PCF]], gene–distance analysis [[TIME_GENES]] and TLS calling and scoring [[TIME_TLS]], with peak memory [[PEAK_MEM]] (Supplementary Table S2). The dominant costs are the repeated `FindMarkers()` calls in the distance profile and the silhouette computation in the parameter search. Users can skip the former by requesting a subset of metrics, and can bound the latter with `max_sample`.

## 4 Discussion

stGradient provides a reproducible workflow for asking boundary-centric questions of single-cell spatial transcriptomics data. Its main design choice is to reduce every analysis to a common coordinate, the signed distance to the nearest tumor boundary, and to make the construction of that boundary explicit and tunable. This differs from domain-segmentation and neighborhood-enrichment approaches, which describe tissue organization without reference to a particular interface. It also complements them: domains or niches identified by tools such as BANKSY or Giotto can be overlaid on stGradient's distance coordinate, and stGradient polygons can be exported as sf objects to any GIS-aware tool.

The breast cancer example shows how combining the modules produces a more specific hypothesis than any one of them alone. Distance profiles indicated that the compositional interface is sharp, within about 30 µm. Cell-type trajectories showed that it is layered, with a Macrophages_1-rich cuff close to the tumor and T cells enriched farther out. PCFs showed that CD8^+^ T cells are specifically associated with these macrophages beyond what shared stromal habitat explains. Gene gradients and TLS calls confirmed that the distal stroma is immunologically active. Taken together, the data are consistent with peritumoral macrophages holding CD8^+^ T cells in the stroma. This hypothesis can be tested directly with replicate sections, cell-state analysis within macrophages and T cells, and functional studies.

Several limitations deserve emphasis, and the package documentation discusses each in detail. First, **cells are not independent replicates.** *P* values from per-cell Wilcoxon tests and regressions are extremely small because *n* is large and neighboring cells are spatially autocorrelated. They should be treated as rankings. Biological claims require replication with the patient or section as the experimental unit (Squair *et al.*, 2021). Second, **the boundary inherits the quality of the malignant-cell annotation**, and it is a two-dimensional cut through a three-dimensional structure, so in-plane distances are upper bounds on true distances. Third, **the PCF null and observation window matter.** By default, the window is the bounding box of the section. This includes tumor interiors and empty background where no extratumoral cells can occur, which inflates *g*(*r*) for all pairs. The CSR null also treats tissue as homogeneous. For this reason, we emphasize relative comparisons (e.g. CD8–macrophage against CD8–stromal) and recommend supplying a window equal to the tissue minus the tumor polygons, as well as a random-labelling null for interaction claims. A future release will make both the default. [[ADJUST IF IMPLEMENTED BEFORE SUBMISSION]] Fourth, **composition and cell state are confounded** in pseudobulk metrics and in gene–distance trends. Segmentation spill-over at the interface produces artifactual gradients in tumor-specific genes, as seen for KRT7, EPCAM and FOXA1. Within-cell-type analyses address the former, and caution about canonical tumor markers addresses the latter. Fifth, **cumulative distance shells smooth the profiles**, which dilutes localized effects. Annular bins would localize effects more sharply and are planned. Finally, the TLS signature and the structural call are not independent, because both are dominated by B- and T-cell markers. The mature-like flag is lenient, since a single LAMP3^+^ DC suffices. Stricter maturity calls should use follicular dendritic cell or germinal-center markers where the panel allows.

Planned extensions include native SpatialExperiment input, support for additional imaging platforms (CosMx, MERSCOPE, Visium HD), multi-sample summaries that treat the patient as the unit of analysis, annular distance bins, and a random-labelling PCF null with tissue-masked windows. stGradient is open source, and we welcome contributions.

## Funding

[[FUNDING]]

*Conflict of interest:* [[COI — e.g. "None declared."]]

## Data availability

The Xenium FFPE human breast cancer dataset analyzed in this article is publicly available from 10x Genomics ([[DATASET_URL]]) and is described in Janesick *et al.* (2023). Code to reproduce all figures (`vignette_figures.Rmd`) is included in the stGradient repository.

## Acknowledgements

[[ACKNOWLEDGEMENTS]]

---

## Table 1. Main exported functions of stGradient, grouped by module

| Module | Function | Purpose |
|---|---|---|
| Setup | `load_xenium_seurat()` | Load object, enforce unique IDs, derive `is_cancer`, cache coordinates |
| Region definition | `launch_stGradient()` | Shiny/plotly app for drawing tumor polygons |
| | `optimize_dbscan_params()` | Grid search over `eps`/`minPts`; silhouette, coverage, purity, F~β~ |
| | `add_dbscan_polygons()`, `get_cancer_regions()` | DBSCAN + concave hull polygons per nest; store/retrieve |
| | `add_polygon_membership()`, `plot_dbscan_polygons()` | Per-cell ROI membership; visual QC |
| Distance | `calculate_boundary_distances()`, `compute_distance_from_boundary()` | Signed distance to nearest boundary |
| Distance profiles | `calculate_distance_profile()` | *D*~tx~, *D*~comp~ and DE strength over distance thresholds |
| | `find_critical_distance()` | Maximum and plateau (90%) distance |
| | `calculate_transcriptomic_distance()`, `calculate_composition_difference()`, `calculate_de_strength()` | Individual metrics for arbitrary cell groups |
| Spatial association | `calculate_pcf_outside_summary()` | Cross-type PCF at distance *r*, extratumoral cells only |
| | `plot_pcf_outside_heatmap()`, `plot_pcf_outside_pair_envelopes()`, `plot_pcf_pair_gg()` | PCF heatmap; full curves with CSR envelopes |
| Gene–distance | `find_distance_associated_genes()`, `plot_distance_genes()` | Per-gene OLS + Spearman vs. distance; volcano and trends |
| TLS | `call_tls()`, `get_tls_summary()` | DBSCAN-based TLS calling; per-aggregate composition and maturity |
| | `add_tls_score()`, `summarize_tls_score()`, `tls_signature_genes()` | 22-gene module score; in- vs. out-of-TLS comparison |
| | `plot_tls_map()`, `plot_tls_composition()`, `plot_tls_score()`, `plot_tls_score_comparison()` | TLS visualization |
| Visualization | `qc_plot_xenium()`, `plot_roi_*()`, `plot_all_metrics()`, `plot_celltype_heatmap()`, `plot_celltype_proportions()` | QC, compartment composition and distance-profile plots |

---

## Figure legends

**Fig. 1. Tumor region definition.** (**A**) The interactive Shiny selector. Cells are colored by malignant annotation (red), and a user-drawn polygon (blue) is returned to R as an sf object. (**B**) Malignant-cell annotations (top, red) and DBSCAN/concave-hull polygons (bottom, blue) for the Xenium breast cancer section; coordinates in µm. (**C**) Schematic of the coverage–purity trade-off. A tight neighborhood radius yields pure polygons that miss scattered tumor cells, whereas a loose radius captures all tumor cells but also stroma. (**D**) Grid search over `eps` and `minPts`, showing coverage, purity and F~0.5~ for combinations with < 15% noise; the selected combination is outlined in red.

**Fig. 2. Distance-resolved composition and expression at the tumor–stroma interface.** (**A**) Cell counts by type inside and outside tumor polygons. (**B**) Composition difference, DE strength and transcriptomic distance as a function of the outer distance threshold *d*. Blue: whole-tumor reference (`"all"`); orange: boundary-band reference (`"distance"`). (**C**) Per-cell-type contribution to composition difference at each threshold. (**D**) Proportion of selected cell types inside the tumor (orange; negative distances) and in cumulative extratumoral shells (blue; positive distances). Black lines are smoothed fits with 95% bands. Macrophages_1 peak within about 50 µm of the boundary, whereas CD4^+^ and CD8^+^ T cells keep increasing with distance.

**Fig. 3. CD8^+^ T cells co-localize with peritumoral macrophages.** (**A**) Cross-type pair correlation *g*(*r*) at *r* = [[PCF_R]] µm among extratumoral cells; values > 1 indicate co-localization. The default bounding-box window inflates absolute values, so the matrix should be read relatively (see Discussion). (**B**) Full PCF for CD8^+^ T cells vs. Macrophages_1 (black) with the pointwise CSR envelope from [[NSIM]] simulations (grey) and the theoretical value 1 (red dashed). (**C**) The same analysis for CD8^+^ T cells vs. stromal cells, a shared-habitat reference pair. Note the different x-axis ranges in B and C.

**Fig. 4. Genes associated with distance from the tumor boundary in peritumoral stroma.** (**A**) Volcano plot of OLS slope (expression change per µm) against −log~10~ BH-adjusted *P*. Triangles mark 88 genes with adjusted *P* below machine precision, drawn at the top. (**B**) Mean expression of top genes in ten quantile-based distance bins.

**Fig. 5. Tertiary lymphoid structures.** (**A**) Nineteen called TLS on the tissue; cells colored by TLS identifier, B cells in red. (**B**) Cellular composition of each TLS. (**C**) Per-cell TLS signature score (top) and with called-TLS perimeters outlined (bottom). (**D**) TLS signature score in cells inside vs. outside called TLS (two-sided Wilcoxon test; cells are not independent replicates).

---

## References

Baddeley,A. *et al.* (2015) *Spatial Point Patterns: Methodology and Applications with R*. Chapman and Hall/CRC, Boca Raton, FL.

Cabrita,R. *et al.* (2020) Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature*, **577**, 561–565.

Cable,D.M. *et al.* (2022) Robust decomposition of cell type mixtures in spatial transcriptomics. *Nat. Biotechnol.*, **40**, 517–526.

Chen,D.S. and Mellman,I. (2017) Elements of cancer immunity and the cancer–immune set point. *Nature*, **541**, 321–330.

Dries,R. *et al.* (2021) Giotto: a toolbox for integrative analysis and visualization of spatial expression data. *Genome Biol.*, **22**, 78.

Ester,M. *et al.* (1996) A density-based algorithm for discovering clusters in large spatial databases with noise. In: *Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96)*. AAAI Press, pp. 226–231.

Hahsler,M. *et al.* (2019) dbscan: fast density-based clustering with R. *J. Stat. Softw.*, **91**, 1–30.

Hao,Y. *et al.* (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat. Biotechnol.*, **42**, 293–304.

Helmink,B.A. *et al.* (2020) B cells and tertiary lymphoid structures promote immunotherapy response. *Nature*, **577**, 549–555.

Hunter,M.V. *et al.* (2021) Spatially resolved transcriptomics reveals the architecture of the tumor-microenvironment interface. *Nat. Commun.*, **12**, 6278.

Janesick,A. *et al.* (2023) High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. *Nat. Commun.*, **14**, 8353.

Keren,L. *et al.* (2018) A structured tumor-immune microenvironment in triple negative breast cancer revealed by multiplexed ion beam imaging. *Cell*, **174**, 1373–1387.

Mariathasan,S. *et al.* (2018) TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. *Nature*, **554**, 544–548.

Moses,L. and Pachter,L. (2022) Museum of spatial transcriptomics. *Nat. Methods*, **19**, 534–546.

Palla,G. *et al.* (2022) Squidpy: a scalable framework for spatial omics analysis. *Nat. Methods*, **19**, 171–178.

Park,J.-S. and Oh,S.-J. (2012) A new concave hull algorithm and concaveness measure for n-dimensional datasets. *J. Inf. Sci. Eng.*, **28**, 587–600.

Pebesma,E. (2018) Simple features for R: standardized support for spatial vector data. *R J.*, **10**, 439–446.

Peranzoni,E. *et al.* (2018) Macrophages impede CD8 T cells from reaching tumor cells and limit the efficacy of anti–PD-1 treatment. *Proc. Natl. Acad. Sci. USA*, **115**, E4041–E4050.

Petitprez,F. *et al.* (2020) B cells are associated with survival and immunotherapy response in sarcoma. *Nature*, **577**, 556–560.

Righelli,D. *et al.* (2022) SpatialExperiment: infrastructure for spatially-resolved transcriptomics data in R using Bioconductor. *Bioinformatics*, **38**, 3128–3131.

Rousseeuw,P.J. (1987) Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *J. Comput. Appl. Math.*, **20**, 53–65.

Sautès-Fridman,C. *et al.* (2019) Tertiary lymphoid structures in the era of cancer immunotherapy. *Nat. Rev. Cancer*, **19**, 307–325.

Sievert,C. (2020) *Interactive Web-Based Data Visualization with R, plotly, and shiny*. Chapman and Hall/CRC, Boca Raton, FL.

Singhal,V. *et al.* (2024) BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. *Nat. Genet.*, **56**, 431–441.

Squair,J.W. *et al.* (2021) Confronting false discoveries in single-cell differential expression. *Nat. Commun.*, **12**, 5692.

Wickham,H. (2016) *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag, New York.
