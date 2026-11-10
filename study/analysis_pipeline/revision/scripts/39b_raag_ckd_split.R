# Renal age-associated genes (RAAG) split by concordance with CKD: are the age genes that also change in CKD the
# reproducible part of the kidney age signal? Same gene-level test and matched null as the benchmark.
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
src <- readLines(file.path(REV, "scripts", "38_gtex_multitissue.R")); eval(parse(text = src[1:grep("^## ---- GTEx tissues", src) - 1]))
set.seed(20261106)
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); Z <- I$Z; grp <- sub("_.*", "", names(Z))
pooled_z <- function(use) { g <- Reduce(union, lapply(Z[use], names)); M <- sapply(Z[use], function(z) z[g]); rownames(M) <- g
  k <- rowSums(is.finite(M)); s <- rowSums(M, na.rm = TRUE) / sqrt(pmax(k, 1)); s[k < max(2, floor(length(use) / 2))] <- NA; s }
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda"))
raag <- unique(lib$gene[lib$sig_id == "A08_Perco2019_RAAG634"]); ckd <- lib[lib$sig_id == "A09_Perco2019_RAAG_CKD", ]; ckd <- ckd[!duplicated(ckd$gene), ]
S2 <- list("RAAG, CKD-concordant" = setNames(rep(0, nrow(ckd)), ckd$gene), "RAAG, not CKD-concordant" = setNames(rep(0, length(setdiff(raag, ckd$gene))), setdiff(raag, ckd$gene)))
K <- do.call(rbind, lapply(names(L <- list("All eligible datasets" = c("G6", "TCGA", "CPTAC"), "Without tumour-adjacent tissue" = c("G6", "TCGA"))), function(a) {
  z <- pooled_z(names(Z)[grp %in% L[[a]]]); z <- z[is.finite(z)]; b <- kab_bins[intersect(names(z), names(kab_bins))]; z <- z[names(b)]
  do.call(rbind, lapply(names(S2), function(k) { g <- intersect(names(S2[[k]]), names(z)); pool <- split(names(b), b)
    nul <- replicate(2000, mean(abs(z[vapply(g, function(x) { p <- pool[[b[[x]]]]; p[sample.int(length(p), 1)] }, "")])))
    data.frame(analysis = a, set = k, n = length(g), ER = mean(abs(z[g])) / mean(nul), lo = mean(abs(z[g])) / quantile(nul, .975), hi = mean(abs(z[g])) / quantile(nul, .025), p = (1 + sum(nul >= mean(abs(z[g])))) / 2001) })) }))
print(K, digits = 3, row.names = FALSE); write.csv(K, file.path(O, "REVISION_raag_ckd_split.csv"), row.names = FALSE)
