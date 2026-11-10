# Extension 1: M-gene for 9 additional C-class gene sets, identical code/null to script 21; then combined C-class summary
suppressMessages({library(metafor)}); REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R")); O <- file.path(REV, "benchmark", "results"); set.seed(20261101)
src <- readLines(file.path(REV, "scripts", "21_benchmark.R")); eval(parse(text=src[1:grep("^Z <- lapply", src)-1]))
I <- readRDS(path.expand("~/kidney_aging_data/benchmark_intermediates.rds")); PZ <- I$PZ
EXTF <- Sys.getenv("EXTF", "signature_library_extension1.csv"); ext <- read.csv(file.path(REV, "benchmark", "signatures", EXTF)); sx <- split(ext, ext$sig_id)
mu <- rowMeans(D$CPTAC$e); sdv <- apply(D$CPTAC$e, 1, sd); bin <- setNames(paste(cut(rank(mu), 10, labels=FALSE), cut(rank(sdv), 10, labels=FALSE)), names(mu))
pz <- PZ$none; ok <- names(bin)[names(bin) %in% names(pz)[is.finite(pz)]]; pool_by_bin <- split(ok, bin[ok])
MG <- do.call(rbind, lapply(names(sx), function(id){ g <- intersect(intersect(unique(sx[[id]]$gene), names(pz)[is.finite(pz)]), names(bin)); if (length(g) < 2) return(NULL)
  obs <- mean(abs(pz[g])); nul <- replicate(2000, mean(abs(pz[vapply(g, function(x) { pl <- pool_by_bin[[bin[x]]]; pl[sample.int(length(pl), 1)] }, "")])))
  data.frame(sig_id=id, n_measured=length(g), mean_abs_z=obs, null_mean=mean(nul), enrichment=obs/mean(nul), emp_p=(1+sum(nul >= obs))/2001, pct_nominal=100*mean(abs(pz[g]) > 1.96)) }))
B <- read.csv(file.path(O, "BENCHMARK_summary.csv")); allp <- c(B$emp_p[is.finite(B$emp_p)], MG$emp_p); MG$fdr_all43 <- p.adjust(allp, "BH")[length(allp) - nrow(MG) + seq_len(nrow(MG))]
B$fdr_all43 <- NA; B$fdr_all43[is.finite(B$emp_p)] <- p.adjust(allp, "BH")[seq_len(sum(is.finite(B$emp_p)))]
lab <- unique(ext[, c("sig_id","class","label","pmid")]); MG <- merge(lab, MG, by="sig_id")
write.csv(MG, file.path(O, "M_gene_extension1.csv"), row.names=FALSE); print(format(MG, digits=2), row.names=FALSE)
Call <- rbind(B[B$class=="C", c("sig_id","label","pmid","n_measured","enrichment","emp_p","fdr_all43","pct_nominal")], MG[, c("sig_id","label","pmid","n_measured","enrichment","emp_p","fdr_all43","pct_nominal")])
cat("\nAll C-class (n =", nrow(Call), "): median enrichment", round(median(Call$enrichment), 3), "| FDR<0.05:", sum(Call$fdr_all43 < .05), "| enrichment<1:", sum(Call$enrichment < 1), "\n")
write.csv(Call[order(-Call$enrichment), ], file.path(O, "C_class_all25.csv"), row.names=FALSE); write.csv(B, file.path(O, "BENCHMARK_summary_fdr43.csv"), row.names=FALSE)
