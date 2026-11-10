# Post hoc revision analyses (2026-11-01), not part of the locked plan:
#  (1) co-primary gene-level benchmark without tumour-adjacent tissue (G6 + TCGA), with the expression-matched null
#      re-derived from the evaluation datasets themselves instead of CPTAC; the locked analysis is re-run alongside;
#  (2) power of the gene-level test as a function of gene-set size and true enrichment ratio, simulated from the
#      observed gene-level statistics and the expression bins occupied by class C genes.
suppressMessages({library(data.table)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "10_cptac_helpers.R")); O <- file.path(REV, "benchmark", "results")
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); Z <- I$Z
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda")); load(file.path(REV, "kidneyAgeBench/data/kab_library.rda"))
grp <- setNames(sub("_.*", "", names(Z)), names(Z))
## ---- expression matrices of the evaluation datasets (same processing as script 21)
E <- list(); cb <- readRDS(path.expand("~/kidney_aging_data/harmonized/harmonized_bundle_corrected.rds"))
for (g in names(cb$harm)) E[[paste0("G6_", g)]] <- cb$harm[[g]]$expr
x <- fread(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_tpm.tsv"), data.table = FALSE); tm <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_meta.csv"))
sx <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_case_sex.csv")); tm$sex <- sx$sex[match(tm$patient, sx$patient)]; tm$age <- ifelse(is.finite(tm$age_days), tm$age_days / 365.25, tm$age_at_index)
tm <- tm[order(tm$sample), ]; tm <- tm[!duplicated(tm$patient) & is.finite(tm$age) & !is.na(tm$sex), ]
e <- log2(as.matrix(x[, tm$file_id]) + 1); o <- order(rowMeans(e), decreasing = TRUE); e <- e[o, ]; gg <- x$gene_name[o]; e <- e[!duplicated(gg), ]; rownames(e) <- gg[!duplicated(gg)]; E$TCGA <- e[rowMeans(e > 0) >= 0.2, ]
V <- path.expand("~/kidney_aging_data/validation/cptac"); mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
tr <- cptac_rd(V, "bcm-ccrcc-transcriptomics"); tr <- cbind(gene = tr[[1]], tr[, grep("_A$", names(tr))]); R <- cptac_sym(tr, e2s); colnames(R) <- sub("_A$", "", colnames(R)); R <- R[rowMeans(R > 0) >= 0.2, ]
cs <- intersect(colnames(R), C$case[is.finite(C$age) & C$sex %in% c("Male", "Female")]); E$CPTAC <- R[, cs]
rm(x, e, R, tr, cb); invisible(gc())
## ---- helpers
min_k <- function(n) max(2, floor(n / 2))
pooled_z <- function(use) { g <- Reduce(union, lapply(Z[use], names)); M <- sapply(Z[use], function(z) z[g]); rownames(M) <- g
  k <- rowSums(is.finite(M)); s <- rowSums(M, na.rm = TRUE) / sqrt(pmax(k, 1)); s[k < min_k(length(use))] <- NA; s }
# expression bins from the evaluation data: within-dataset percentile of mean and SD, averaged across datasets, then deciles
eval_bins <- function(use) { P <- lapply(E[use], function(m) { mu <- rowMeans(m); sdv <- apply(m, 1, sd); cbind(mu = rank(mu) / length(mu), sd = rank(sdv) / length(sdv)) })
  g <- Reduce(union, lapply(P, rownames)); mu <- sapply(P, function(p) p[, "mu"][g]); sdv <- sapply(P, function(p) p[, "sd"][g]); rownames(mu) <- rownames(sdv) <- g
  k <- rowSums(is.finite(mu)); keep <- k >= min_k(length(use)); m1 <- rowMeans(mu[keep, ], na.rm = TRUE); s1 <- rowMeans(sdv[keep, ], na.rm = TRUE)
  setNames(paste(cut(rank(m1), 10, labels = FALSE), cut(rank(s1), 10, labels = FALSE)), names(m1)) }
leak <- function(id) if (grepl("^A1[12]_", id)) "G6" else if (grepl("^A0[56]_Rowland", id)) "TCGA" else character()
null_means <- function(g, pz, bins, pool, n) { M <- vapply(g, function(x) { p <- pool[[bins[[x]]]]; abs(pz[p[sample.int(length(p), n, TRUE)]]) }, numeric(n)); if (is.null(dim(M))) mean(M) else rowMeans(M) }
bench <- function(sets, groups, bins, nnull = 2000) do.call(rbind, lapply(sets, function(s) { id <- s$sig_id[1]; use <- names(Z)[grp %in% setdiff(groups, leak(id))]
  pz <- pooled_z(use); ok <- intersect(names(pz)[is.finite(pz)], names(bins)); g <- intersect(unique(s$gene), ok); if (length(g) < 2) return(NULL)
  pool <- split(ok, bins[ok]); obs <- mean(abs(pz[g])); nul <- null_means(g, pz, bins, pool, nnull)
  data.frame(sig_id = id, class = s$class[1], n = length(g), enrichment = obs / mean(nul), emp_p = (1 + sum(nul >= obs)) / (nnull + 1)) }))
## ---- (1) four configurations: data {with, without CPTAC} x null bins {CPTAC (locked), evaluation data}
sets <- split(kab_library, kab_library$sig_id)
cfg <- list(primary_locked = list(g = c("G6", "TCGA", "CPTAC"), b = kab_bins), primary_evalbins = list(g = c("G6", "TCGA", "CPTAC"), b = NULL),
            noCPTAC_cptacbins = list(g = c("G6", "TCGA"), b = kab_bins), noCPTAC_evalbins = list(g = c("G6", "TCGA"), b = NULL))
for (k in names(cfg)) if (is.null(cfg[[k]]$b)) cfg[[k]]$b <- eval_bins(names(E)[grp[names(E)] %in% cfg[[k]]$g])
set.seed(20261101)
RES <- do.call(rbind, lapply(names(cfg), function(k) { r <- bench(sets, cfg[[k]]$g, cfg[[k]]$b); r$fdr <- p.adjust(r$emp_p, "BH"); r$config <- k; r }))
RES$status <- ifelse(RES$fdr < .05, "FDR < 0.05", ifelse(RES$enrichment >= 1, "n.s., ER >= 1", "ER < 1"))
write.csv(RES, file.path(O, "REVISION_gene_level_four_configs.csv"), row.names = FALSE)
SUM <- do.call(rbind, lapply(split(RES, list(RES$config, RES$class)), function(d) data.frame(config = d$config[1], class = d$class[1], n_sets = nrow(d),
  median_ER = median(d$enrichment), q1 = quantile(d$enrichment, .25), q3 = quantile(d$enrichment, .75), n_fdr05 = sum(d$fdr < .05), n_ER_below1 = sum(d$enrichment < 1))))
wil <- sapply(names(cfg), function(k) { d <- RES[RES$config == k, ]; wilcox.test(d$enrichment[d$class == "A"], d$enrichment[d$class == "C"])$p.value })
SUM$wilcox_A_vs_C_p <- wil[SUM$config]; write.csv(SUM, file.path(O, "REVISION_class_summary_four_configs.csv"), row.names = FALSE); print(SUM, row.names = FALSE, digits = 3)
cat("Spearman of ER, locked vs co-primary:", round(cor(RES$enrichment[RES$config == "primary_locked"], RES$enrichment[RES$config == "noCPTAC_evalbins"][match(RES$sig_id[RES$config == "primary_locked"], RES$sig_id[RES$config == "noCPTAC_evalbins"])], method = "spearman", use = "complete.obs"), 2), "\n")
## pooled gene pools against the matched null, per configuration (class C genes carry no leakage)
lib <- kab_library; ga <- unique(lib$gene[lib$sig_id == "B03_GenAge"]); ca <- unique(lib$gene[lib$sig_id == "B04_CellAge"])
pools <- list("Class A, external sets" = unique(lib$gene[lib$class == "A" & !grepl("^A1[12]_|^A0[56]_", lib$sig_id)]), GenAge = ga, CellAge = ca, "Class C" = unique(lib$gene[lib$class == "C"]))
POOL <- do.call(rbind, lapply(names(cfg), function(k) { use <- names(Z)[grp %in% cfg[[k]]$g]; pz <- pooled_z(use); b <- cfg[[k]]$b; ok <- intersect(names(pz)[is.finite(pz)], names(b)); pool <- split(ok, b[ok])
  do.call(rbind, lapply(names(pools), function(p) { g <- intersect(pools[[p]], ok); obs <- mean(abs(pz[g])); nul <- null_means(g, pz, b, pool, 2000)
    data.frame(config = k, pool = p, n_genes = length(g), naive_ratio = obs / mean(abs(pz[ok])), matched_ratio = obs / mean(nul), lo = obs / quantile(nul, .975), hi = obs / quantile(nul, .025), emp_p = (1 + sum(nul >= obs)) / 2001) })) }))
write.csv(POOL, file.path(O, "REVISION_gene_pools_four_configs.csv"), row.names = FALSE); print(POOL, row.names = FALSE, digits = 3)
## ---- (2) power simulation
# A true enrichment ratio is induced by sampling member genes within their expression bins with weights exp(theta * |z|);
# theta is calibrated so that the expected ratio over the class C bin mixture equals the target. theta = 0 gives the null.
power_sim <- function(k_cfg, ks = c(2, 3, 4, 5, 6, 8, 10, 13, 20, 30, 50, 100), ers = c(1, 1.1, 1.2, 1.3, 1.5, 2), reps = 1000, nnull = 1000) {
  use <- names(Z)[grp %in% cfg[[k_cfg]]$g]; pz <- pooled_z(use); b <- cfg[[k_cfg]]$b; ok <- intersect(names(pz)[is.finite(pz)], names(b))
  az <- lapply(split(ok, b[ok]), function(g) abs(pz[g])); tb <- table(b[intersect(pools[["Class C"]], ok)]); w <- setNames(as.numeric(tb) / sum(tb), names(tb))
  mb <- function(th) sapply(names(w), function(n) { v <- az[[n]]; sum(v * exp(th * v)) / sum(exp(th * v)) }); base <- sum(w * mb(0))
  th <- sapply(ers, function(r) if (r == 1) 0 else uniroot(function(t) sum(w * mb(t)) / base - r, c(0, 10))$root)
  cw <- lapply(th, function(t) lapply(az, function(v) cumsum(exp(t * v))))
  out <- list()
  for (k in ks) { hit <- matrix(NA, reps, length(ers)); erl <- matrix(NA, reps, length(ers))
    for (r in seq_len(reps)) { bb <- sample(names(w), k, TRUE, prob = w)
      nul <- rowMeans(vapply(bb, function(n) { v <- az[[n]]; v[sample.int(length(v), nnull, TRUE)] }, numeric(nnull)))
      for (j in seq_along(ers)) { obs <- mean(vapply(bb, function(n) { cc <- cw[[j]][[n]]; az[[n]][findInterval(runif(1) * cc[length(cc)], cc) + 1] }, 0))
        hit[r, j] <- (1 + sum(nul >= obs)) / (nnull + 1) < .05; erl[r, j] <- obs / mean(nul) } }
    out[[length(out) + 1]] <- data.frame(config = k_cfg, k = k, true_ER = ers, theta = th, power = colMeans(hit), mean_realised_ER = colMeans(erl), median_realised_ER = apply(erl, 2, median)) }
  do.call(rbind, out) }
set.seed(20261102)
PW <- rbind(power_sim("noCPTAC_evalbins"), power_sim("primary_locked"))
write.csv(PW, file.path(O, "REVISION_power_by_setsize.csv"), row.names = FALSE); print(PW, row.names = FALSE, digits = 3)
## smallest true ER detectable with 80% power at each class C set size (log-linear interpolation over the simulated grid)
mde <- function(d, k) { p <- sapply(sort(unique(d$true_ER)), function(r) { s <- d[d$true_ER == r, ]; approx(log(s$k), s$power, log(k), rule = 2)$y }); er <- sort(unique(d$true_ER))
  if (max(p) < .8) return(Inf); if (p[1] >= .8) return(er[1]); i <- which(p >= .8)[1]; er[i - 1] + (er[i] - er[i - 1]) * (.8 - p[i - 1]) / (p[i] - p[i - 1]) }
MDE <- do.call(rbind, lapply(c("noCPTAC_evalbins", "primary_locked"), function(k) { d <- RES[RES$config == k & RES$class == "C", ]; pw <- PW[PW$config == k, ]
  data.frame(config = k, sig_id = d$sig_id, n = d$n, enrichment = d$enrichment, fdr = d$fdr, mde80 = sapply(d$n, function(n) mde(pw, n)),
             power_ER1.3 = sapply(d$n, function(n) { s <- pw[pw$true_ER == 1.3, ]; approx(log(s$k), s$power, log(n), rule = 2)$y })) }))
write.csv(MDE, file.path(O, "REVISION_classC_detectable_effect.csv"), row.names = FALSE)
for (k in unique(MDE$config)) { d <- MDE[MDE$config == k, ]; cat(k, ": class C sets with >= 80% power for ER 1.3:", sum(d$power_ER1.3 >= .8), "of", nrow(d), "; median power at ER 1.3:", round(median(d$power_ER1.3), 2), "\n") }
