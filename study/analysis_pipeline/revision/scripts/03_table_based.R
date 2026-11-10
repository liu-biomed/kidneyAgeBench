# Reanalyses that only need the released per-donor / per-sample output tables
suppressMessages({library(metafor); library(pROC)})
T <- "~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/tables"
O <- "~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030/results/tables"
rd <- function(f) read.csv(file.path(path.expand(T), f), stringsAsFactors=FALSE)
wr <- function(x, f) write.csv(x, file.path(path.expand(O), f), row.names=FALSE)
set.seed(20261030)

## 1. Immune-nucleus fraction vs age in the snRNA-seq normal donors (Reviewer 1: 7-donor analysis)
sn <- rd("sn_donor_composition_normal.csv"); sn$is_cortex <- sn$is_cortex %in% c("True", TRUE)
sp_perm <- function(x, y, B=1e5){ r0 <- cor(x, y, method="s"); rp <- replicate(B, cor(x, sample(y), method="s")); c(rho=r0, p_perm=mean(abs(rp) >= abs(r0)-1e-12)) }
boot_ci <- function(x, y, B=1e4){ b <- replicate(B, { i <- sample(seq_along(x), replace=TRUE); suppressWarnings(cor(x[i], y[i], method="s")) }); quantile(b, c(.025,.975), na.rm=TRUE) }
cx <- sn[sn$is_cortex, ]; an <- sn
rows <- list()
for (nm in c("cortex_dominant","all_normal")) {
  d <- if (nm=="cortex_dominant") cx else an
  pp <- sp_perm(d$age, d$imm_frac); ci <- boot_ci(d$age, d$imm_frac)
  rows[[nm]] <- data.frame(subset=nm, n=nrow(d), rho=pp["rho"], perm_p=pp["p_perm"], boot_lo=ci[1], boot_hi=ci[2],
                           age_range=paste(range(d$age), collapse="-"), n_female=sum(d$sex=="female"))
}
imm <- do.call(rbind, rows)
lodo <- do.call(rbind, lapply(seq_len(nrow(cx)), function(i) data.frame(dropped=cx$donor_id[i], dropped_age=cx$age[i],
          rho=cor(cx$age[-i], cx$imm_frac[-i], method="s"), p=cor.test(cx$age[-i], cx$imm_frac[-i], method="s", exact=TRUE)$p.value)))
wr(imm, "R1_sn_immune_fraction_summary.csv"); wr(lodo, "R1_sn_immune_fraction_leave_one_donor_out.csv")
cat("== snRNA immune fraction\n"); print(imm, row.names=FALSE); print(lodo, row.names=FALSE)

## 2. Compartment distribution of age-down and age-up genes (not only PT)
ga <- rd("sn_gene_celltype_assignment.csv")
comp <- as.data.frame.matrix(table(ga$direction, ga$compartment)); comp$direction <- rownames(comp)
pt_share <- aggregate(pt_cpm/(max_cpm) ~ direction, ga, median)
wr(comp, "R1_sn_compartment_by_direction.csv"); cat("== compartments\n"); print(comp)

## 3. Within-cell-type pseudobulk restricted to normal donors and disease-adjusted
pbc <- rd("deconv_intrinsic_pseudobulk_class.csv"); pbs <- rd("deconv_intrinsic_pseudobulk_subclass.csv")
pbs$class_s <- pbs$subclass; pb <- rbind(pbc[,c("donor_id","class_s","age","disease","sex","n","up","dn","sig")],
                                        pbs[pbs$subclass=="PT", c("donor_id","class_s","age","disease","sex","n","up","dn","sig")])
pb_res <- do.call(rbind, lapply(split(pb, pb$class_s), function(d) do.call(rbind, lapply(c("dn","up","sig"), function(a){
  dn <- d[d$disease=="normal",]; f <- lm(scale(d[[a]]) ~ scale(d$age) + d$disease)
  data.frame(cell=d$class_s[1], axis=a, n_all=nrow(d), rho_all=cor(d[[a]], d$age, method="s"),
    p_all=cor.test(d[[a]], d$age, method="s", exact=FALSE)$p.value,
    beta_age_disease_adj=coef(f)[2], p_age_disease_adj=summary(f)$coef[2,4],
    n_normal=nrow(dn), rho_normal=cor(dn[[a]], dn$age, method="s"), p_normal=cor.test(dn[[a]], dn$age, method="s", exact=FALSE)$p.value) }))))
wr(pb_res, "R1_pseudobulk_normal_only_and_disease_adjusted.csv"); cat("== pseudobulk\n"); print(format(pb_res, digits=2), row.names=FALSE)

## 4. Deconvolution: fraction-age pooled across all 7 cohorts vs RNA-seq-only vs microarray-only
fa <- rd("deconv_fraction_age_percohort.csv"); rnaseq <- c("GSE93069","GSE217427")
fa_res <- do.call(rbind, lapply(split(fa, paste(fa$method, fa$celltype)), function(d) do.call(rbind, lapply(c("all","rnaseq","array"), function(s){
  dd <- switch(s, all=d, rnaseq=d[d$cohort %in% rnaseq,], array=d[!d$cohort %in% rnaseq,]); dd <- dd[is.finite(dd$zi) & is.finite(dd$sei),]; if (nrow(dd)==0) return(NULL)
  f <- rma(yi=dd$zi, sei=dd$sei, method=if (nrow(dd)>2) "REML" else "FE")
  data.frame(method=d$method[1], celltype=d$celltype[1], subset=s, k=nrow(dd), n=sum(dd$n), rho=tanh(f$b), lo=tanh(f$ci.lb), hi=tanh(f$ci.ub), p=f$pval) }))))
wr(fa_res, "R1_deconv_fraction_age_by_platform.csv")
cat("== deconvolution PT/immune by platform\n"); print(format(fa_res[fa_res$celltype %in% c("PT","IMM","epithelial","immune"),], digits=2), row.names=FALSE)

## 5. Clinical cohorts: AUC with DeLong 95% CI, arms, and diagnosis-level contrasts (no imposed ordering)
cs <- rd("clin_persample_scores.csv")
cs$grp <- ifelse(cs$cohort=="GSE175759", ifelse(cs$diagnosis=="Control","Reference","Disease"),
          ifelse(cs$cohort=="GSE104954", ifelse(grepl("Reference", cs$group),"Reference","Disease"), cs$group))
auc_rows <- list()
for (co in unique(cs$cohort)) for (a in c("combined_score","up_score","down_score")) {
  d <- cs[cs$cohort==co,]; r <- roc(d$grp=="Disease", d[[a]], direction="<", quiet=TRUE); ci <- ci.auc(r, method="delong")
  auc_rows[[length(auc_rows)+1]] <- data.frame(cohort=co, axis=a, n_ref=sum(d$grp=="Reference"), n_dis=sum(d$grp=="Disease"),
     AUC=as.numeric(auc(r)), lo=ci[1], hi=ci[3]) }
auc_t <- do.call(rbind, auc_rows)
# combined vs up-arm alone (DeLong paired)
for (co in unique(cs$cohort)) { d <- cs[cs$cohort==co,]
  r1 <- roc(d$grp=="Disease", d$combined_score, direction="<", quiet=TRUE); r2 <- roc(d$grp=="Disease", d$up_score, direction="<", quiet=TRUE)
  auc_t$p_combined_vs_up[auc_t$cohort==co] <- roc.test(r1, r2, method="delong")$p.value }
# random-effects pooled AUC (logit scale) to show heterogeneity
ac <- auc_t[auc_t$axis=="combined_score",]; lg <- qlogis(ac$AUC); se <- (qlogis(pmin(ac$hi,.999)) - qlogis(pmax(ac$lo,.001)))/(2*1.96)
pa <- rma(yi=lg, sei=se, method="REML"); auc_pool <- data.frame(pooled_AUC=plogis(pa$b), lo=plogis(pa$ci.lb), hi=plogis(pa$ci.ub), I2=pa$I2, Q_p=pa$QEp, pi_lo=plogis(predict(pa)$pi.lb), pi_hi=plogis(predict(pa)$pi.ub))
wr(auc_t, "R1_clinical_AUC_DeLong.csv"); wr(auc_pool, "R1_clinical_AUC_pooled.csv")
cat("== AUC\n"); print(format(auc_t, digits=3), row.names=FALSE); print(auc_pool)
# GSE104954 / GSE175759: each diagnosis vs reference (Hedges g, Wilcoxon), no ordering assumption
dx <- list()
for (co in c("GSE104954","GSE175759")) { d <- cs[cs$cohort==co,]; ref <- d$combined_score[d$grp=="Reference"]
  for (g in setdiff(unique(d$diagnosis[d$grp=="Disease"]), NA)) { x <- d$combined_score[d$grp=="Disease" & d$diagnosis==g]; if (length(x) < 3) next
    sp <- sqrt(((length(x)-1)*var(x)+(length(ref)-1)*var(ref))/(length(x)+length(ref)-2))
    dx[[length(dx)+1]] <- data.frame(cohort=co, diagnosis=g, n=length(x), n_ref=length(ref), hedges_g=(mean(x)-mean(ref))/sp*(1-3/(4*(length(x)+length(ref))-9)),
      wilcox_p=wilcox.test(x, ref, exact=FALSE)$p.value) } }
dx <- do.call(rbind, dx); dx$fdr <- p.adjust(dx$wilcox_p, "BH")
kw <- sapply(c("GSE104954","GSE175759"), function(co){ d <- cs[cs$cohort==co & !is.na(cs$diagnosis),]; kruskal.test(d$combined_score, factor(d$diagnosis))$p.value })
wr(dx, "R1_clinical_diagnosis_vs_reference.csv"); cat("== diagnosis contrasts\n"); print(format(dx, digits=2), row.names=FALSE); print(kw)
# GSE66494 arm distributions
g66 <- cs[cs$cohort=="GSE66494",]; print(aggregate(cbind(up_score, down_score, combined_score) ~ grp, g66, function(v) round(c(mean=mean(v), sd=sd(v)),2)))
