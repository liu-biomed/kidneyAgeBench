#' Gene-level age statistics in normal human kidney (unadjusted for composition)
#'
#' List of named numeric vectors (Fisher z of partial Spearman correlation with age) for six
#' GEO donor cohorts (GSE36895, GSE60807, GSE147451 donor-level, GSE66272, GSE93069, GSE167093),
#' TCGA KIRC/KIRP/KICH solid tissue normal (n = 129), CPTAC ccRCC normal adjacent RNA (n = 75)
#' and GTEx v8 kidney cortex (n = 85). Attribute "group" gives G6/TCGA/CPTAC/GTEx.
#' @format A list of 9 named numeric vectors.
"kab_gene_z"
#' Composition-adjusted gene-level age statistics in normal human kidney
#'
#' List of named numeric vectors (Fisher z of partial Spearman correlation with age, adjusted
#' for sex/project and five marker-based tissue composition scores: immune, stromal, endothelial,
#' distal tubular and proximal tubular content). Attribute "group" gives G6/TCGA/CPTAC/GTEx.
#' @format A list of 9 named numeric vectors.
"kab_gene_z_comp"
#' Expression-matched bins (CPTAC NAT RNA mean x SD deciles) used for the locked primary analysis.
#' @format Named character vector.
"kab_bins"
#' Expression-matched bins (G6 + TCGA mean x SD deciles) used for the co-primary analysis.
#' @format Named character vector.
"kab_bins_coprimary"
#' Simulated statistical power grid by gene-set size and true enrichment ratio.
#' @format data.frame with columns config, k, true_ER, power.
"kab_power_grid"
#' Library of 48 published kidney-ageing, ageing and senescence gene sets (classes A, B, C).
#' @format data.frame with columns sig_id, class, label, pmid, gene, dir.
"kab_library"
