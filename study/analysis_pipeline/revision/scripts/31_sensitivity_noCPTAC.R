# Sensitivity: gene-level benchmark restricted to G6 + TCGA (CPTAC removed); same code, null and leakage rules as scripts 21/23/24
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
O <- file.path(REV, "benchmark", "results"); set.seed(20261103)
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); Z <- I$Z
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda")); load(file.path(REV, "kidneyAgeBench/data/kab_library.rda"))
grp <- sub("_.*", "", names(Z))
pz_fun <- function(excl) { use <- names(Z)[grp %in% setdiff(c("G6","TCGA"), excl)]; g <- Reduce(union, lapply(Z[use], names))
  M <- sapply(Z[use], function(z) z[g]); rownames(M) <- g; k <- rowSums(is.finite(M))
  s <- rowSums(M, na.rm=TRUE)/sqrt(pmax(k,1)); s[k < max(2, floor(length(use)/2))] <- NA; s }
leak <- function(id) if (grepl("^A1[12]_", id)) "G6" else if (grepl("^A0[56]_Rowland", id)) "TCGA" else character()
R <- do.call(rbind, lapply(split(kab_library, kab_library$sig_id), function(s){ id <- s$sig_id[1]; pz <- pz_fun(leak(id))
  ok <- intersect(names(pz)[is.finite(pz)], names(kab_bins)); g <- intersect(unique(s$gene), ok); if (length(g) < 2) return(NULL)
  pool <- split(ok, kab_bins[ok]); obs <- mean(abs(pz[g]))
  nul <- replicate(2000, mean(abs(pz[vapply(g, function(x) { p <- pool[[kab_bins[[x]]]]; p[sample.int(length(p), 1)] }, "")])))
  data.frame(sig_id=id, class=s$class[1], n=length(g), enrichment=obs/mean(nul), emp_p=(1+sum(nul >= obs))/2001) }))
R$fdr <- p.adjust(R$emp_p, "BH"); write.csv(R, file.path(O, "SENS_gene_level_G6_TCGA_only.csv"), row.names=FALSE)
for (k in c("A","B","C")) { d <- R[R$class==k,]; cat(sprintf("class %s: n=%d median enrichment=%.2f FDR<0.05=%d below1=%d\n", k, nrow(d), median(d$enrichment), sum(d$fdr<.05), sum(d$enrichment<1))) }
cat("GenAge", round(R$enrichment[R$sig_id=="B03_GenAge"],2), "CellAge", round(R$enrichment[R$sig_id=="B04_CellAge"],2), "\n")
cat("note: A11/A12 omitted in gene-level G6+TCGA analysis (only 1 non-overlapping dataset, <2 required); A05/A06 evaluated in G6 (6 cohorts)\n")
