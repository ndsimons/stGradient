# Suppress R CMD check notes for NSE/data-masked symbols used in dplyr/ggplot2.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".ct", ".data", "celltype", "cluster", "composition_difference",
    "contribution", "ct_grp", "de_strength", "distance", "inside_mode",
    "inside_n", "inside_prop", "inside_polygon", "is_cancer", "lbl", "n",
    "n_type1", "n_type2", "outside_prop", "prop", "proportion", "region",
    "roi_status", "total_n", "transcriptomic_distance", "value", "x", "y",
    "cellType", "g_value", "type1", "type2", "pCR", "patient_id",
    "niche", "niche_label", "wilcox_p", "perm_p", "p_value", "p_adj",
    "p_adj_lm", "mean_pCR", "mean_non_pCR", "perm_padj", "group",
    "mean_prop", "sem", "sig_label", "estimate", "std_error", "coef",
    "direction", "gene", "bin_midpoint", "mean_expr",
    "n_B", "n_T", "n_DC", "n_total", "frac_B", "is_TLS", "has_matureDC",
    "in_TLS", "TLS_id", "eps", "minPts", "metric", "neg_log_p", "capped"
  ))
}
