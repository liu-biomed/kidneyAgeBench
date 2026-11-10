# Reviewer 1 point 3: fully nested leave-one-cohort-out validation. In each fold, the genome-wide meta-analysis,
# FDR control, sign-consistency, coverage and LODO filters, direction assignment and scoring are rebuilt from the
# training cohorts only; the resulting score is then evaluated against age in the held-out cohort (and in GTEx).
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
b <- load_bundle(); harm <- b$harm; gt <- b$gtex; P <- readRDS(MP)
Z <- P$Z; SE <- P$SE; cohorts <- colnames(Z)
eval_sig <- function(up, dn, h, adjust_sex){
  if (length(up) + length(dn) == 0) return(c(rho=NA, p=NA, prho=NA, n_genes=0))
  s <- score_z(h$expr, up, dn); ok <- is.finite(s$score); a <- h$meta$age
  sp <- spear(s$score[ok], a[ok]); pr <- NA
  if (adjust_sex && !all(is.na(h$meta$sex))) pr <- cor(residuals(lm(rank(s$score[ok]) ~ h$meta$sex[ok])), residuals(lm(rank(a[ok]) ~ h$meta$sex[ok])))
  c(rho=unname(sp[1]), p=unname(sp[2]), prho=pr, n_genes=s$n_up[1] + s$n_dn[1]) }
folds <- list(); sigs <- list()
for (h in cohorts) {
  tr <- setdiff(cohorts, h); Zt <- Z[, tr]; St <- SE[, tr]
  res <- meta_genes(Zt, St, PAR$mink_loco)
  rob <- select_robust(res, Zt, St, fdr=0.05, sc=0.7, kmin=PAR$kmin_loco)   # coverage rule scaled to the number of training cohorts
  rule <- "prespecified"
  if (nrow(rob) < 10) {                                                 # fallback: fixed-size signature
    cand <- res[res$sign_consist >= 0.7 & res$k >= PAR$kmin_loco, ]; cand <- cand[order(cand$pval), ][1:nrow(P$rob), ]
    cand$direction <- ifelse(cand$z > 0, "up", "down"); rob <- cand; rule <- "top94_fallback" }
  up <- rob$gene[rob$direction=="up"]; dn <- rob$gene[rob$direction=="down"]
  # always also evaluate a fixed-size top-94 signature for comparability across folds
  c94 <- res[res$sign_consist >= 0.7 & res$k >= PAR$kmin_loco, ]; c94 <- c94[order(c94$pval), ][1:nrow(P$rob), ]
  up94 <- c94$gene[c94$z > 0]; dn94 <- c94$gene[c94$z < 0]
  hh <- harm[[h]]; sx <- !(h %in% NO_SEX)
  e1 <- eval_sig(up, dn, hh, sx); e2 <- eval_sig(up94, dn94, hh, sx)
  e3 <- eval_sig(up, dn, gt, TRUE)
  folds[[h]] <- data.frame(heldout=h, n=ncol(hh$expr), train_fdr05=sum(res$fdr < 0.05), sig_rule=rule, sig_size=nrow(rob),
    overlap_with_full94=length(intersect(rob$gene, P$rob$gene)), heldout_genes_measured=e1["n_genes"],
    rho=e1["rho"], p=e1["p"], partial_rho_sex=e1["prho"], rho_top94=e2["rho"], p_top94=e2["p"],
    gtex_rho_fold_signature=e3["rho"], gtex_p=e3["p"])
  sigs[[h]] <- list(up=up, dn=dn)
  cat(h, "done:", nrow(rob), "genes, held-out rho =", round(e1["rho"], 3), "\n")
}
F <- do.call(rbind, folds); print(format(F, digits=3), row.names=FALSE); wr(F, "R1_nested_LOCO_folds.csv")
# In-sample (non-nested) rho for the full 94-gene signature, to quantify selection-induced optimism
ins <- sapply(cohorts, function(h) eval_sig(P$rob$gene[P$rob$direction=="up"], P$rob$gene[P$rob$direction=="down"], harm[[h]], FALSE)["rho"])
pooled <- rbind(cbind(analysis="nested LOCO (prespecified rule)", pool_rho(F$rho, F$n)),
                cbind(analysis="nested LOCO (fixed top-94)", pool_rho(F$rho_top94, F$n)),
                cbind(analysis="in-sample full 94-gene signature", pool_rho(ins, F$n)))
print(pooled); wr(pooled, "R1_nested_LOCO_pooled.csv")
saveRDS(list(folds=F, sigs=sigs, insample=ins), file.path(D, "harmonized", paste0("loco_", KSET, ".rds")))
