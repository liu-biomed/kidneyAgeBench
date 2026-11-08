#' Pooled gene-level age statistic across normal-kidney datasets
#'
#' Stouffer combination of per-dataset Fisher-z age statistics (partial Spearman with
#' sex, and project for TCGA; optionally adjusted for 5 marker-based tissue-composition scores).
#' Genes measured in fewer than half of the used datasets (minimum 2) are set to NA.
#' @param groups Dataset groups to use: any of "G6", "TCGA", "CPTAC", "GTEx".
#' @param exclude Groups to drop (e.g. those used to derive the signature under test).
#' @param composition_adjusted Logical; if TRUE, uses composition-adjusted gene-level age z statistics.
#' @return Named numeric vector of pooled z statistics (positive = higher with age).
#' @export
kab_pooled_z <- function(groups = c("G6", "TCGA", "CPTAC"), exclude = character(), composition_adjusted = FALSE) {
  Z <- if (isTRUE(composition_adjusted)) kidneyAgeBench::kab_gene_z_comp else kidneyAgeBench::kab_gene_z
  use <- names(Z)[attr(Z, "group") %in% setdiff(groups, exclude)]
  if (!length(use)) stop("no datasets left after exclusion")
  g <- Reduce(union, lapply(Z[use], names))
  M <- sapply(Z[use], function(z) z[g]); if (is.null(dim(M))) M <- matrix(M, ncol = 1, dimnames = list(g, use))
  rownames(M) <- g; k <- rowSums(is.finite(M))
  s <- rowSums(M, na.rm = TRUE) / sqrt(pmax(k, 1)); s[k < max(2, floor(length(use) / 2))] <- NA
  s
}

#' Benchmark one gene set against age in normal human kidney (Implements Box 1 Checks 1-5)
#'
#' @param genes Character vector of HGNC symbols.
#' @param direction Optional numeric vector (same length as `genes`): +1 expected to rise
#'   with age, -1 expected to fall, 0 or NA unsigned.
#' @param exclude Dataset groups to exclude (leakage control), e.g. "TCGA" for TCGA-derived sets.
#' @param config Either "locked" (G6 + TCGA + CPTAC with CPTAC bins) or "coprimary" (G6 + TCGA with evaluation-data bins).
#' @param n_null Number of expression-matched random gene sets.
#' @param seed Random seed.
#' @return A one-row data.frame: genes measured, mean |z|, null mean, enrichment ratio and empirical P
#'   (both unadjusted and composition-adjusted), estimated power for ER = 1.3 and minimum detectable ER
#'   at 80% power (`mde80`), percentage of members with |z| > 1.96, and directional concordance with
#'   both expression-matched signed-null P (`D_excess_z`, `p_dir_matched`) and binomial P (`p_dir`).
#' @export
kab_benchmark <- function(genes, direction = NULL, exclude = character(),
                          config = c("locked", "coprimary"), n_null = 2000, seed = 1) {
  config <- match.arg(config)
  grps <- if (config == "locked") c("G6", "TCGA", "CPTAC") else c("G6", "TCGA")
  bins <- if (config == "locked") kidneyAgeBench::kab_bins else kidneyAgeBench::kab_bins_coprimary
  pz <- kab_pooled_z(groups = grps, exclude = exclude, composition_adjusted = FALSE)
  pzc <- kab_pooled_z(groups = grps, exclude = exclude, composition_adjusted = TRUE)
  if (is.null(direction)) direction <- rep(0, length(genes))
  d <- data.frame(gene = toupper(genes), dir = direction)[!duplicated(toupper(genes)), ]
  ok <- names(pz)[is.finite(pz) & is.finite(pzc)]; ok <- intersect(ok, names(bins)); d <- d[d$gene %in% ok, ]
  if (nrow(d) < 2) stop("fewer than 2 genes measured")
  pool <- split(ok, bins[ok]); set.seed(seed)
  idx <- vapply(d$gene, function(x) { p <- pool[[bins[[x]]]]; p[sample.int(length(p), n_null, TRUE)] }, character(n_null))
  Zn <- matrix(pz[idx], nrow = n_null)
  Znc <- matrix(pzc[idx], nrow = n_null)
  obs <- mean(abs(pz[d$gene])); nul <- rowMeans(abs(Zn))
  obsc <- mean(abs(pzc[d$gene])); nulc <- rowMeans(abs(Znc))
  s <- !is.na(d$dir) & d$dir != 0; hit <- sign(pz[d$gene[s]]) == sign(d$dir[s])
  if (any(s)) {
    D_obs <- mean(d$dir[s] * pz[d$gene[s]])
    D_nul <- as.numeric(Zn[, s, drop = FALSE] %*% d$dir[s]) / sum(s)
    D_excess_z <- (D_obs - mean(D_nul)) / stats::sd(D_nul)
    p_dir_matched <- (1 + sum(D_nul >= D_obs)) / (n_null + 1)
    D_obsc <- mean(d$dir[s] * pzc[d$gene[s]])
    D_nulc <- as.numeric(Znc[, s, drop = FALSE] %*% d$dir[s]) / sum(s)
    D_excess_z_comp <- (D_obsc - mean(D_nulc)) / stats::sd(D_nulc)
    p_dir_matched_comp <- (1 + sum(D_nulc >= D_obsc)) / (n_null + 1)
  } else {
    D_excess_z <- NA_real_; p_dir_matched <- NA_real_
    D_excess_z_comp <- NA_real_; p_dir_matched_comp <- NA_real_
  }
  pw_cfg <- if (config == "locked") "primary_locked" else "noCPTAC_evalbins"
  pw <- kidneyAgeBench::kab_power_grid[kidneyAgeBench::kab_power_grid$config == pw_cfg, ]
  pw13 <- stats::approx(log(pw$k[pw$true_ER == 1.3]), pw$power[pw$true_ER == 1.3], log(nrow(d)), rule = 2)$y
  ers <- sort(unique(pw$true_ER))
  p_grid <- sapply(ers, function(r) stats::approx(log(pw$k[pw$true_ER == r]), pw$power[pw$true_ER == r], log(nrow(d)), rule = 2)$y)
  mde80 <- if (max(p_grid) < 0.8) Inf else if (p_grid[1] >= 0.8) ers[1] else {
    i <- which(p_grid >= 0.8)[1]; ers[i - 1] + (ers[i] - ers[i - 1]) * (0.8 - p_grid[i - 1]) / (p_grid[i] - p_grid[i - 1])
  }
  data.frame(n_measured = nrow(d), mean_abs_z = obs, null_mean = mean(nul), enrichment = obs / mean(nul),
             emp_p = (1 + sum(nul >= obs)) / (n_null + 1),
             enrichment_comp = obsc / mean(nulc), emp_p_comp = (1 + sum(nulc >= obsc)) / (n_null + 1),
             power_ER1.3 = pw13, mde80 = mde80,
             pct_nominal = 100 * mean(abs(pz[d$gene]) > 1.96),
             n_signed = sum(s), pct_dir_correct = if (any(s)) 100 * mean(hit) else NA,
             D_excess_z = D_excess_z, p_dir_matched = p_dir_matched,
             D_excess_z_comp = D_excess_z_comp, p_dir_matched_comp = p_dir_matched_comp,
             p_dir = if (any(s)) stats::binom.test(sum(hit), sum(s), alternative = "greater")$p.value else NA)
}

#' Re-run the benchmark for every gene set in the bundled library
#' @param n_null Random sets per signature.
#' @param config Either "locked" or "coprimary".
#' @return data.frame with one row per gene set and BH-adjusted FDR.
#' @export
kab_benchmark_library <- function(n_null = 2000, config = c("locked", "coprimary")) {
  config <- match.arg(config)
  L <- kidneyAgeBench::kab_library; out <- lapply(split(L, L$sig_id), function(s) {
    ex <- if (grepl("^A1[12]_", s$sig_id[1])) "G6" else if (grepl("^A0[56]_Rowland", s$sig_id[1])) "TCGA" else character()
    r <- tryCatch(kab_benchmark(s$gene, s$dir, exclude = ex, config = config, n_null = n_null), error = function(e) NULL)
    if (is.null(r)) return(NULL); cbind(s[1, c("sig_id", "class", "label", "pmid")], r) })
  R <- do.call(rbind, out); R$fdr <- stats::p.adjust(R$emp_p, "BH"); R$fdr_comp <- stats::p.adjust(R$emp_p_comp, "BH")
  R[order(R$class, -R$enrichment), ]
}
