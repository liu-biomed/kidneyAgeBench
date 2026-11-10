# Addendum 2 D3 (same statistic and permutation scheme as 11_decoupling.R; vectorised implementation for speed).
# Per-gene Spearman(protein, mRNA) within younger vs older halves; complete pairs only, computed by ranks column-wise.
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "10_cptac_helpers.R"))
set.seed(20261030)
V <- path.expand("~/kidney_aging_data/validation/cptac"); mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
P <- cptac_sym(cptac_rd(V, "bcm-ccrcc-proteomics.*Normal"), e2s); P <- P[rowMeans(is.finite(P)) >= 0.5, ]
tr <- cptac_rd(V, "bcm-ccrcc-transcriptomics"); tr <- cbind(gene=tr[[1]], tr[, grep("_A$", names(tr))]); R <- cptac_sym(tr, e2s); colnames(R) <- sub("_A$", "", colnames(R)); R <- R[rowMeans(R > 0) >= 0.2, ]
cs <- Reduce(intersect, list(colnames(P), colnames(R), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")])); m <- C[match(cs, C$case), ]
gb <- intersect(rownames(P), rownames(R)); Pm <- P[gb, cs]; Rm <- R[gb, cs]
# restrict to genes complete in both layers (Spearman on complete data makes vectorisation exact)
ok <- rowSums(is.finite(Pm)) == ncol(Pm) & rowSums(is.finite(Rm)) == ncol(Rm); Pm <- Pm[ok,]; Rm <- Rm[ok,]; cat("complete genes:", nrow(Pm), "| patients:", ncol(Pm), "\n")
rowspear <- function(A, B){ A <- t(apply(A, 1, rank)); B <- t(apply(B, 1, rank)); A <- A - rowMeans(A); B <- B - rowMeans(B); rowSums(A*B)/sqrt(rowSums(A^2)*rowSums(B^2)) }
stat <- function(yg){ ry <- rowspear(Pm[, yg], Rm[, yg]); ro <- rowspear(Pm[, !yg], Rm[, !yg]); cl <- function(x) atanh(pmin(pmax(x, -.999), .999))
  c(median_diff=median(cl(ro) - cl(ry), na.rm=TRUE), young=median(ry, na.rm=TRUE), old=median(ro, na.rm=TRUE)) }
young <- m$age <= median(m$age); obs <- stat(young); nul <- replicate(2000, stat(sample(young))["median_diff"])
d3 <- data.frame(n_young=sum(young), n_old=sum(!young), age_split=median(m$age), genes=nrow(Pm), young_median_r=obs["young"], old_median_r=obs["old"],
                 median_fisher_diff_old_minus_young=obs["median_diff"], perm_p_two_sided=(1+sum(abs(nul) >= abs(obs["median_diff"])))/2001)
print(d3); write.csv(d3, file.path(O, "D3_decoupling.csv"), row.names=FALSE)
