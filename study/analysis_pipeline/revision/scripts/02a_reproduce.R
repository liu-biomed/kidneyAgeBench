# Reproduce the original genome-wide meta-analysis and 94-gene robust signature from rebuilt cohorts
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
b <- load_bundle(); harm <- b$harm
for (g in names(harm)) cat(sprintf("%-10s genes=%5d n=%3d sex=%s\n", g, nrow(harm[[g]]$expr), ncol(harm[[g]]$expr), !(g %in% NO_SEX)))
eff <- effects_all(harm); M <- zse(eff)
t0 <- Sys.time(); res <- meta_genes(M$Z, M$SE, PAR$mink); cat("meta genes:", nrow(res), "time", round(difftime(Sys.time(), t0, units="mins"),1), "min\n")
rob <- select_robust(res, M$Z, M$SE, kmin=PAR$kmin)
orig <- read.csv(file.path(REV, "..", "tables", "meta_analysis_full_results.csv"))
orob <- read.csv(file.path(REV, "..", "tables", "robust_aging_signature.csv"))
mm <- merge(res, orig, by="gene")
repro <- data.frame(metric=c("genes meta-analysed (orig 18025)", "genes FDR<0.05 (orig 136)", "robust genes (orig 94)",
  "robust overlap with original 94", "cor pooled r vs original", "max |pooled r diff|", "genes in original not in rebuild"),
  value=c(nrow(res), sum(res$fdr<0.05), nrow(rob), length(intersect(rob$gene, orob$gene)), cor(mm$r, mm$pooled_r),
          max(abs(mm$r-mm$pooled_r)), length(setdiff(orig$gene, res$gene))))
print(repro); wr(repro, "R0_reproduction_check.csv")
cat("missing from rebuild:", paste(setdiff(orob$gene, rob$gene), collapse=","), " | new:", paste(setdiff(rob$gene, orob$gene), collapse=","), "\n")
saveRDS(list(eff=eff, Z=M$Z, SE=M$SE, res=res, rob=rob), MP)
