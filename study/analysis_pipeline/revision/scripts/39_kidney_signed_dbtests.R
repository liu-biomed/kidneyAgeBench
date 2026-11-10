# Post hoc (2026-11-01): the signed gene-level test of script 38, applied in the kidney evaluation data to every signed
# general resource and to the CellAge inducer/inhibitor annotation, so kidney and GTEx tissues are compared on one scale.
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
src <- readLines(file.path(REV, "scripts", "38_gtex_multitissue.R"))
eval(parse(text = src[1:grep("^## ---- GTEx tissues", src) - 1]))   # gene sets and helper functions only
set.seed(20261105)
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); Z <- I$Z; grp <- sub("_.*", "", names(Z))
pooled_z <- function(use) { g <- Reduce(union, lapply(Z[use], names)); M <- sapply(Z[use], function(z) z[g]); rownames(M) <- g
  k <- rowSums(is.finite(M)); s <- rowSums(M, na.rm = TRUE) / sqrt(pmax(k, 1)); s[k < max(2, floor(length(use) / 2))] <- NA; s }
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda"))
sets <- c("CellAge inducers (+) / inhibitors (-)", "de Magalhães 2009", "Palmer 2021", "SenMayo", "Fridman 2008", "GenAge", "CellAge")
K <- do.call(rbind, lapply(names(L <- list("Locked primary datasets" = c("G6", "TCGA", "CPTAC"), "Co-primary datasets" = c("G6", "TCGA"))), function(a) {
  z <- pooled_z(names(Z)[grp %in% L[[a]]]); z <- z[is.finite(z)]; b <- kab_bins[intersect(names(z), names(kab_bins))]; z <- z[names(b)]
  do.call(rbind, lapply(sets, function(k) cbind(analysis = a, set = k, test_set(SETS[[k]], z, b)))) }))
K$fdr_signed <- ave(K$p_signed, K$analysis, FUN = function(p) p.adjust(p, "BH"))
write.csv(K, file.path(O, "REVISION_kidney_signed_dbtests.csv"), row.names = FALSE)
print(K[, c("analysis", "set", "n", "ER", "p_unsigned", "D_excess_z", "p_signed", "fdr_signed", "pct_concordant")], digits = 3, row.names = FALSE)
