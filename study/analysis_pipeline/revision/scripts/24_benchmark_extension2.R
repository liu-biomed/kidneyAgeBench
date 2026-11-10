# Extension 1: M-gene for 9 additional C-class gene sets, identical code/null to script 21; then combined C-class summary
suppressMessages({library(metafor)}); REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "benchmark", "results"); set.seed(20261102)
src <- readLines(file.path(REV, "scripts", "21_benchmark.R")); eval(parse(text=src[1:grep("^Z <- lapply", src)-1]))
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); PZ <- I$PZ
EXTF <- Sys.getenv("EXTF", "signature_library_extension1.csv"); ext <- read.csv(file.path(REV, "benchmark", "signatures", EXTF)); sx <- split(ext, ext$sig_id)
mu <- rowMeans(D$CPTAC$e); sdv <- apply(D$CPTAC$e, 1, sd); bin <- setNames(paste(cut(rank(mu), 10, labels=FALSE), cut(rank(sdv), 10, labels=FALSE)), names(mu))
pz <- PZ$none; ok <- names(bin)[names(bin) %in% names(pz)[is.finite(pz)]]; pool_by_bin <- split(ok, bin[ok])
MG <- do.call(rbind, lapply(names(sx), function(id){ g <- intersect(intersect(unique(sx[[id]]$gene), names(pz)[is.finite(pz)]), names(bin)); if (length(g) < 2) return(NULL)
  obs <- mean(abs(pz[g])); nul <- replicate(2000, mean(abs(pz[vapply(g, function(x) { pl <- pool_by_bin[[bin[x]]]; pl[sample.int(length(pl), 1)] }, "")])))
  data.frame(sig_id=id, n_measured=length(g), mean_abs_z=obs, null_mean=mean(nul), enrichment=obs/mean(nul), emp_p=(1+sum(nul >= obs))/2001, pct_nominal=100*mean(abs(pz[g]) > 1.96)) }))
lab <- unique(ext[, c("sig_id","class","label","pmid")]); MG <- merge(lab, MG, by="sig_id"); write.csv(MG, file.path(O, "M_gene_extension2.csv"), row.names=FALSE)
## final combined table over all gene sets (initial 34 + extension 1 + extension 2), BH across all
B <- read.csv(file.path(O, "BENCHMARK_summary.csv")); E1 <- read.csv(file.path(O, "M_gene_extension1.csv"))
cols <- c("sig_id","class","label","pmid","n_measured","enrichment","emp_p","pct_nominal")
A <- rbind(B[is.finite(B$emp_p), cols], E1[, cols], MG[, cols]); A$fdr_all <- p.adjust(A$emp_p, "BH")
A <- merge(A, B[, c("sig_id","pct_dir_correct","pooled_rho","pooled_lo","pooled_hi","pooled_p","pooled_rho_comp","GTEx_rho")], by="sig_id", all.x=TRUE)
A <- A[order(A$class, -A$enrichment), ]; write.csv(A, file.path(O, "BENCHMARK_FINAL_all_sets.csv"), row.names=FALSE)
print(format(MG, digits=2), row.names=FALSE)
for (k in c("A","B","C")) { d <- A[A$class==k,]; cat(sprintf("class %s: n=%d median enrichment=%.2f FDR<0.05=%d enrichment<1=%d\n", k, nrow(d), median(d$enrichment), sum(d$fdr_all < .05), sum(d$enrichment < 1))) }
