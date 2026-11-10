# Post hoc extension (2026-11-01): do ageing-database genes track age outside the kidney?
# In each GTEx v8 tissue (n >= 70 donors), gene-level age statistics are computed exactly as in the kidney benchmark
# (rank-based partial correlation, Fisher z), and gene sets are compared with random genes matched on expression mean
# and SD deciles of the same tissue. Signed tests check the direction claimed by a set, including the CellAge
# inducer/inhibitor annotation (inducers expected to rise with age, inhibitors to fall). The same signed CellAge test is
# also run in the kidney evaluation data (locked primary and co-primary datasets).
suppressMessages({library(data.table)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
O <- file.path(REV, "benchmark", "results"); G <- path.expand("~/kidney_aging_data/gtex"); set.seed(20261104)
load(file.path(REV, "kidneyAgeBench/data/kab_library.rda")); lib <- kab_library
long <- read.csv(file.path(REV, "benchmark/signatures/signature_library_long.csv"))
## ---- gene sets
ca <- read.delim(path.expand("~/kidney_aging_data/published/hagr/cellage3.tsv"))
cl <- long[long$sig_id == "B04_CellAge", c("gene_orig", "gene")]; ca$gene <- cl$gene[match(ca$Gene.symbol, cl$gene_orig)]
ce <- aggregate(Senescence.Effect ~ gene, ca[!is.na(ca$gene) & ca$Senescence.Effect %in% c("Induces", "Inhibits"), ], function(v) paste(sort(unique(v)), collapse = "/"))
ind <- ce$gene[ce$Senescence.Effect == "Induces"]; inh <- ce$gene[ce$Senescence.Effect == "Inhibits"]   # genes annotated both ways are dropped
gs <- function(id) { d <- lib[lib$sig_id == id, ]; d <- d[!duplicated(d$gene), ]; setNames(d$dir, d$gene) }
SETS <- list(GenAge = gs("B03_GenAge"), CellAge = gs("B04_CellAge"),
  "CellAge inducers (+) / inhibitors (-)" = c(setNames(rep(1, length(ind)), ind), setNames(rep(-1, length(inh)), inh)),
  "de Magalhães 2009" = gs("B01_deMagalhaes2009"), "Palmer 2021" = gs("B02_Palmer2021"), SenMayo = gs("B05_SenMayo"), "Fridman 2008" = gs("B06_Fridman2008_up"), "Sood 2015" = gs("B08_Sood2015"))
cg <- unique(lib$gene[lib$class == "C"]); ag <- unique(lib$gene[lib$class == "A" & !grepl("^A1[12]_|^A0[56]_", lib$sig_id)])
SETS[["Class C genes (kidney disease 'ageing' models)"]] <- setNames(rep(0, length(cg)), cg); SETS[["Class A genes (kidney age signatures, external)"]] <- setNames(rep(0, length(ag)), ag)
cat("CellAge inducers:", length(ind), "inhibitors:", length(inh), "dropped (both):", sum(grepl("/", ce$Senescence.Effect)), "\n")
## ---- helpers (same statistics as scripts 21 and 36)
gene_z <- function(e, age, base) { X <- if (is.null(base)) matrix(1, length(age)) else model.matrix(~ ., base); q <- ncol(X) - 1
  Hm <- X %*% solve(crossprod(X)) %*% t(X); a <- rank(age, na.last = "keep"); a <- a - Hm %*% a
  Rk <- t(apply(e, 1, rank, na.last = "keep")); Rk <- Rk - t(Hm %*% t(Rk)); r <- as.numeric((Rk / sqrt(rowSums(Rk^2))) %*% (a / sqrt(sum(a^2))))
  setNames(atanh(pmin(pmax(r, -.9999), .9999)) * sqrt(length(age) - q - 3), rownames(e)) }
bins_of <- function(e) { mu <- rowMeans(e); sdv <- apply(e, 1, sd); setNames(paste(cut(rank(mu), 10, labels = FALSE), cut(rank(sdv), 10, labels = FALSE)), rownames(e)) }
test_set <- function(s, z, b, nnull = 2000) { g <- intersect(names(s), names(z)[is.finite(z)]); if (length(g) < 5) return(NULL)
  pool <- split(names(b), b); dir <- s[g]; sg_mask <- !is.na(dir) & dir != 0; signed <- any(sg_mask)
  idx <- vapply(g, function(x) { p <- pool[[b[[x]]]]; p[sample.int(length(p), nnull, TRUE)] }, character(nnull))
  Zn <- matrix(z[idx], nrow = nnull); obs_u <- mean(abs(z[g])); nul_u <- rowMeans(abs(Zn))
  out <- data.frame(n = length(g), ER = obs_u / mean(nul_u), p_unsigned = (1 + sum(nul_u >= obs_u)) / (nnull + 1), signed = signed)
  if (signed) { D <- mean(dir[sg_mask] * z[g[sg_mask]]); nulD <- as.numeric(Zn[, sg_mask, drop = FALSE] %*% dir[sg_mask]) / sum(sg_mask)
    out$D <- D; out$D_null_mean <- mean(nulD); out$D_null_sd <- sd(nulD); out$D_excess_z <- (D - mean(nulD)) / sd(nulD); out$p_signed <- (1 + sum(nulD >= D)) / (nnull + 1); out$pct_concordant <- 100 * mean(sign(z[g[sg_mask]]) == sign(dir[sg_mask])) } else { out$D <- NA; out$D_null_mean <- NA; out$D_null_sd <- NA; out$D_excess_z <- NA; out$p_signed <- NA; out$pct_concordant <- NA }
  out }
## ---- GTEx tissues
ph <- fread(file.path(G, "SubjectPhenotypesDS.txt"), data.table = FALSE); ph$age <- as.numeric(substr(ph$AGE, 1, 2)) + 4.5
sa <- fread(file.path(G, "SampleAttributesDS.txt"), data.table = FALSE, select = c("SAMPID", "SMTSISCH"))
files <- list.files(file.path(G, "tpm_by_tissue"), pattern = "gct.gz$", full.names = TRUE)
RES <- list(); EDA <- list(); NT <- list()
for (f in files) { tis <- sub("gene_tpm_2017-06-05_v8_(.*)\\.gct\\.gz", "\\1", basename(f))
  x <- fread(cmd = paste("gzcat", shQuote(f)), skip = 2, data.table = FALSE); sym <- x$Description; m <- as.matrix(x[, -(1:3)]); rm(x)
  sub_id <- sub("^([^-]+-[^-]+).*", "\\1", colnames(m)); keep <- !duplicated(sub_id) & sub_id %in% ph$SUBJID; m <- m[, keep]; sub_id <- sub_id[keep]
  if (ncol(m) < 70) { cat(tis, ": n =", ncol(m), "skipped\n"); next }
  e <- log2(m + 1); ok <- rowMeans(m >= 0.1) >= 0.2 & !is.na(sym) & sym != ""; e <- e[ok, ]; s2 <- sym[ok]; o <- order(rowMeans(e), decreasing = TRUE); e <- e[o, ]; s2 <- s2[o]; e <- e[!duplicated(s2), ]; rownames(e) <- s2[!duplicated(s2)]; rm(m)
  p <- ph[match(sub_id, ph$SUBJID), ]; isch <- sa$SMTSISCH[match(colnames(e), sa$SAMPID)]; isch[!is.finite(isch)] <- median(isch, na.rm = TRUE)
  hardy <- factor(ifelse(is.na(p$DTHHRDY), "NA", p$DTHHRDY)); sex <- factor(p$SEX)
  base_full <- data.frame(hardy = hardy, isch = isch); if (nlevels(sex) > 1) base_full$sex <- sex; base_full <- base_full[, sapply(base_full, function(v) length(unique(v)) > 1), drop = FALSE]
  base_sex <- if (nlevels(sex) > 1) data.frame(sex = sex) else NULL
  z_full <- gene_z(e, p$age, base_full); z_sex <- gene_z(e, p$age, base_sex); b <- bins_of(e)
  NT[[tis]] <- data.frame(tissue = tis, n = ncol(e), genes = nrow(e), age_min = min(p$age), age_max = max(p$age))
  EDA[[tis]] <- data.frame(tissue = tis, EDA2R_z = unname(z_full["EDA2R"]), CDKN2A_z = unname(z_full["CDKN2A"]))
  for (adj in c("full", "sex_only")) { z <- if (adj == "full") z_full else z_sex
    RES[[paste(tis, adj)]] <- do.call(rbind, lapply(names(SETS), function(k) { r <- test_set(SETS[[k]], z, b); if (is.null(r)) return(NULL); cbind(tissue = tis, adjustment = adj, set = k, r) })) }
  cat(tis, ": n =", ncol(e), "genes =", nrow(e), "\n"); rm(e); invisible(gc()) }
R <- do.call(rbind, RES); R$fdr_unsigned <- ave(R$p_unsigned, R$set, R$adjustment, FUN = function(p) p.adjust(p, "BH")); R$fdr_signed <- ave(R$p_signed, R$set, R$adjustment, FUN = function(p) p.adjust(p, "BH"))
write.csv(R, file.path(O, "REVISION_gtex_multitissue_sets.csv"), row.names = FALSE)
write.csv(do.call(rbind, NT), file.path(O, "REVISION_gtex_multitissue_tissues.csv"), row.names = FALSE); write.csv(do.call(rbind, EDA), file.path(O, "REVISION_gtex_multitissue_controls.csv"), row.names = FALSE)
S <- do.call(rbind, lapply(split(R, list(R$set, R$adjustment), drop = TRUE), function(d) data.frame(set = d$set[1], adjustment = d$adjustment[1], tissues = nrow(d),
  median_ER = median(d$ER), q1 = quantile(d$ER, .25), q3 = quantile(d$ER, .75), n_enriched_fdr05 = sum(d$fdr_unsigned < .05 & d$ER > 1),
  n_signed_fdr05 = if (all(is.na(d$p_signed))) NA else sum(d$fdr_signed < .05, na.rm = TRUE), median_pct_concordant = median(d$pct_concordant))))
write.csv(S, file.path(O, "REVISION_gtex_multitissue_summary.csv"), row.names = FALSE); print(S[order(S$adjustment, -S$median_ER), ], row.names = FALSE, digits = 3)
E <- do.call(rbind, EDA); cat("EDA2R z > 0 in", sum(E$EDA2R_z > 0, na.rm = TRUE), "of", sum(is.finite(E$EDA2R_z)), "tissues; median z", round(median(E$EDA2R_z, na.rm = TRUE), 2), "\n")
## ---- the signed CellAge test in the kidney evaluation data (locked primary and co-primary datasets)
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); Z <- I$Z; grp <- sub("_.*", "", names(Z))
RV <- read.csv(file.path(O, "REVISION_gene_level_four_configs.csv"))
pooled_z <- function(use) { g <- Reduce(union, lapply(Z[use], names)); M <- sapply(Z[use], function(z) z[g]); rownames(M) <- g
  k <- rowSums(is.finite(M)); s <- rowSums(M, na.rm = TRUE) / sqrt(pmax(k, 1)); s[k < max(2, floor(length(use) / 2))] <- NA; s }
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda"))
K <- do.call(rbind, lapply(list(locked = c("G6", "TCGA", "CPTAC"), coprimary = c("G6", "TCGA")), function(gr) { z <- pooled_z(names(Z)[grp %in% gr]); z <- z[is.finite(z)]
  b <- kab_bins[intersect(names(z), names(kab_bins))]; z <- z[names(b)]
  do.call(rbind, lapply(c("CellAge inducers (+) / inhibitors (-)", "GenAge", "CellAge"), function(k) cbind(set = k, test_set(SETS[[k]], z, b)))) }))
K$analysis <- rep(c("locked datasets, CPTAC bins", "co-primary datasets, CPTAC bins"), each = 3); print(K, row.names = FALSE, digits = 3)
write.csv(K, file.path(O, "REVISION_kidney_cellage_signed.csv"), row.names = FALSE)
