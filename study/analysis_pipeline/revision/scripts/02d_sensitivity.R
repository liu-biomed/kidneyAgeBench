# Reviewer 1 points 7-10: sex, meta-analytic diagnostics, selection thresholds, platform-stratified discovery, probe collapsing
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
suppressMessages({library(fgsea); library(msigdbr)})
set.seed(20261030)
b <- load_bundle(); harm <- b$harm; gt <- b$gtex; P <- readRDS(MP)
Z <- P$Z; SE <- P$SE; res <- P$res; rob <- P$rob; G <- readRDS(file.path(D,"harmonized",paste0("gtex_scores_", KSET, ".rds"))); ge <- G$ge
H <- as.data.frame(msigdbr(species="Homo sapiens", collection="H")); hl <- split(H$gene_symbol, H$gs_name)
gsea_h <- function(rs){ rk <- rs[rs$k >= max(rs$k)-1, ]; st <- sort(setNames(rk$z, rk$gene), decreasing=TRUE)
  as.data.frame(fgseaMultilevel(hl, st, minSize=10, maxSize=500, eps=0))[, c("pathway","NES","padj")] }
gh0 <- gsea_h(res)
gtex_conc <- function(genes, dirz){ g <- intersect(genes, ge$gene); if (!length(g)) return(c(NA,NA))
  c(100*mean(sign(ge$r[match(g, ge$gene)]) == sign(dirz[g])), length(g)) }

## 1. Sex: restrict to the seven cohorts with sex metadata
s7 <- setdiff(colnames(Z), NO_SEX); r7 <- meta_genes(Z[, s7], SE[, s7], PAR$mink_sex); rb7 <- select_robust(r7, Z[, s7], SE[, s7], kmin=PAR$kmin_sex)
mm <- merge(res[, c("gene","z")], r7[, c("gene","z","pval","fdr")], by="gene", suffixes=c("_all","_sex7"))
k94 <- r7[r7$gene %in% rob$gene, ]; k94 <- merge(k94, rob[, c("gene","z")], by="gene", suffixes=c("","_all"))
gh7 <- gsea_h(r7); ghm <- merge(gh0, gh7, by="pathway", suffixes=c("_all","_sex7"))
# Sex-stratified and age-by-sex interaction for the 94-gene score in the seven cohorts
sx <- do.call(rbind, lapply(s7, function(g){ h <- harm[[g]]; sc <- score_z(h$expr, rob$gene[rob$direction=="up"], rob$gene[rob$direction=="down"])$score
  d <- data.frame(sc=as.numeric(scale(sc)), age=h$meta$age/10, sex=factor(h$meta$sex, levels=c("F","M")))
  f <- lm(sc ~ age*sex, d); fm <- if (sum(d$sex=="M") > 5) lm(sc ~ age, d[d$sex=="M",]) else NULL; ff <- if (sum(d$sex=="F") > 5) lm(sc ~ age, d[d$sex=="F",]) else NULL
  gb <- function(m) if (is.null(m)) c(NA,NA) else summary(m)$coef["age", 1:2]
  data.frame(cohort=g, n_F=sum(d$sex=="F"), n_M=sum(d$sex=="M"), int_b=summary(f)$coef["age:sexM",1], int_se=summary(f)$coef["age:sexM",2],
             bF=gb(ff)[1], seF=gb(ff)[2], bM=gb(fm)[1], seM=gb(fm)[2]) }))
pm <- function(b, s){ ok <- is.finite(b) & is.finite(s); f <- rma(yi=b[ok], sei=s[ok], method="REML"); c(est=f$b, lo=f$ci.lb, hi=f$ci.ub, p=f$pval, k=sum(ok)) }
sexsum <- data.frame(
  metric=c("cohorts / samples", "genes FDR<0.05 (sex-available cohorts)", "robust genes (sex-available cohorts, scaled coverage rule)", "overlap with primary signature",
           "genome-wide correlation of pooled z (all vs sex-available cohorts)", "signature genes: same direction", "signature genes: nominal P<0.05",
           "signature genes: FDR<0.05", "signature genes: median |r| ratio (sex-available / all)", "Hallmark NES correlation (all vs sex-available)",
           "Hallmark FDR<0.05 in both analyses"),
  value=c(paste(length(s7), sum(sapply(harm[s7], function(h) ncol(h$expr))), sep=" / "), sum(r7$fdr<0.05), nrow(rb7), length(intersect(rb7$gene, rob$gene)),
          round(cor(mm$z_all, mm$z_sex7), 3), sum(sign(k94$z)==sign(k94$z_all)), sum(k94$pval<0.05), sum(k94$fdr<0.05),
          round(median(abs(tanh(k94$z))/abs(tanh(k94$z_all))), 2), round(cor(ghm$NES_all, ghm$NES_sex7), 3), sum(ghm$padj_all<0.05 & ghm$padj_sex7<0.05)))
sexint <- rbind(data.frame(term="age x sex interaction (SD per decade, M vs F)", t(pm(sx$int_b, sx$int_se))),
                data.frame(term="age effect in females", t(pm(sx$bF, sx$seF))), data.frame(term="age effect in males", t(pm(sx$bM, sx$seM))))
print(sexsum); print(sexint); wr(sexsum, "R1_sex_restricted_summary.csv"); wr(sexint, "R1_sex_interaction.csv"); wr(sx, "R1_sex_interaction_by_cohort.csv")
wr(ghm[order(ghm$padj_all), ], "R1_sex_restricted_hallmark.csv")

## 2. Meta-analytic diagnostics for the 94 genes
dg <- do.call(rbind, mclapply(rob$gene, function(g){ zi <- Z[g,]; sei <- SE[g,]; ok <- is.finite(zi)&is.finite(sei)
  f <- pool1(zi, sei); fk <- pool1(zi, sei, test="knha"); ci <- tryCatch(confint(f)$random, error=function(e) NULL); if (!is.matrix(ci)) ci <- NULL; pr <- predict(f)
  w <- weights(f); names(w) <- names(zi)[ok]; big <- c("GSE147451","GSE167093"); small <- names(which(sapply(harm, function(h) ncol(h$expr)) < 25))
  fb <- pool1(zi[!names(zi) %in% big], sei[!names(sei) %in% big]); fs <- pool1(zi[!names(zi) %in% small], sei[!names(sei) %in% small])
  pz <- 2*pnorm(-abs(zi[ok]/sei[ok]))
  data.frame(gene=g, k=sum(ok), r=tanh(f$b), tau2=f$tau2, tau2_hi=if (is.null(ci)) NA else ci["tau^2","ci.ub"], I2=f$I2, I2_hi=if (is.null(ci)) NA else ci["I^2(%)","ci.ub"],
    pi_lo=tanh(pr$pi.lb), pi_hi=tanh(pr$pi.ub), p_REML=f$pval, p_HKSJ=fk$pval, weight_two_largest=sum(w[names(w) %in% big]),
    r_without_two_largest=tanh(fb$b), p_without_two_largest=fb$pval, r_without_small=tanh(fs$b), p_without_small=fs$pval,
    sign_consist=mean(sign(zi[ok])==sign(as.numeric(f$b))), n_cohorts_nominal=sum(pz < 0.05 & sign(zi[ok])==sign(as.numeric(f$b)))) }, mc.cores=1))
dsum <- data.frame(metric=c("I2 = 0", "upper 95% CI of I2 > 50%", "prediction interval excludes 0", "P_HKSJ < 0.05", "P_HKSJ < 0.001",
   "median weight of two largest cohorts", "same direction & P<0.05 without two largest cohorts", "same direction & P<0.05 without cohorts n<25",
   "sign consistency = 100%", "sign consistency >= 80%", "median cohorts individually nominal (same sign)"),
   value=c(sum(dg$I2 < 1e-6), sum(dg$I2_hi > 50, na.rm=TRUE), sum(sign(dg$pi_lo)==sign(dg$pi_hi)), sum(dg$p_HKSJ < 0.05), sum(dg$p_HKSJ < 0.001),
           round(median(dg$weight_two_largest), 1), sum(sign(dg$r_without_two_largest)==sign(dg$r) & dg$p_without_two_largest < 0.05),
           sum(sign(dg$r_without_small)==sign(dg$r) & dg$p_without_small < 0.05), sum(dg$sign_consist == 1), sum(dg$sign_consist >= 0.8), median(dg$n_cohorts_nominal)))
print(dsum); wr(dsum, "R1_meta_diagnostics_summary.csv"); wr(dg, "R1_meta_diagnostics_by_gene.csv")

## 3. Threshold sensitivity
im_sets <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_INTERFERON_ALPHA_RESPONSE")
me_sets <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_ADIPOGENESIS")
imm <- unique(unlist(hl[im_sets])); met <- unique(unlist(hl[me_sets])); U <- res$gene[res$k >= PAR$kmin - 1]
hyp <- function(x, S) { x <- intersect(x, U); S <- intersect(S, U); c(ov=length(intersect(x, S)), p=phyper(length(intersect(x,S))-1, length(S), length(U)-length(S), length(x), lower.tail=FALSE)) }
grid <- expand.grid(fdr=c(0.01, 0.05, 0.10), sc=c(0.6, 0.7, 0.8, 1.0), kmin=PAR$kgrid, lodo=c(TRUE, FALSE))
dz <- setNames(res$z, res$gene)
th <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i){ g <- grid[i,]
  s <- select_robust(res, Z, SE, fdr=g$fdr, sc=g$sc, kmin=g$kmin, lodo=g$lodo); u <- s$gene[s$direction=="up"]; d <- s$gene[s$direction=="down"]
  gc <- gtex_conc(s$gene, dz); hi <- hyp(u, imm); hm <- hyp(d, met)
  data.frame(g, n=nrow(s), n_up=length(u), n_down=length(d), jaccard_vs_primary=length(intersect(s$gene, rob$gene))/length(union(s$gene, rob$gene)),
             gtex_concordance_pct=gc[1], up_in_immune_hallmarks=hi[1], p_up_immune=hi[2], down_in_metabolic_hallmarks=hm[1], p_down_metabolic=hm[2]) }))
print(format(th, digits=2), row.names=FALSE); wr(th, "R1_threshold_sensitivity.csv")

## 4. Platform-stratified discovery and cross-platform testing (no gene reuse between selection and test)
arr <- setdiff(colnames(Z), RNASEQ)
ra <- meta_genes(Z[, arr], SE[, arr], PAR$mink_arr); sa <- select_robust(ra, Z[, arr], SE[, arr], kmin=PAR$kmin_arr)
rq <- if (length(RNASEQ) > 1) meta_genes(Z[, RNASEQ], SE[, RNASEQ], 2) else { z <- Z[, RNASEQ]; s <- SE[, RNASEQ]; ok <- is.finite(z)
  data.frame(gene=rownames(Z)[ok], z=z[ok], se=s[ok], pval=2*pnorm(-abs(z[ok]/s[ok])), k=1, sign_consist=1, fdr=p.adjust(2*pnorm(-abs(z[ok]/s[ok])), "BH")) }
test_in <- function(sel, rt){ m <- merge(sel[, c("gene","z")], rt[, c("gene","z","pval")], by="gene", suffixes=c("_sel","_test"))
  c(n_sel=nrow(sel), n_tested=nrow(m), concordance=100*mean(sign(m$z_sel)==sign(m$z_test)), r=cor(m$z_sel, m$z_test),
    binom_p=binom.test(sum(sign(m$z_sel)==sign(m$z_test)), nrow(m))$p.value, nominal_same_dir=sum(sign(m$z_sel)==sign(m$z_test) & m$pval<0.05)) }
sq <- rq[rq$k==length(RNASEQ) & rq$sign_consist==1, ]; sq <- sq[order(sq$pval), ][1:100, ]
plat <- rbind(cbind(direction=paste0("microarray-selected (", length(arr), " cohorts) -> RNA-seq (", length(RNASEQ), ")"), t(test_in(sa, rq))),
              cbind(direction=paste0("RNA-seq top-100 (", length(RNASEQ), ") -> microarray (", length(arr), " cohorts)"), t(test_in(sq, ra))))
# score-level: array-derived signature applied to the RNA-seq cohorts and GTEx
sc_lvl <- do.call(rbind, lapply(c(RNASEQ, "GTEx"), function(g){ h <- if (g=="GTEx") gt else harm[[g]]
  s <- score_z(h$expr, sa$gene[sa$direction=="up"], sa$gene[sa$direction=="down"]); sp <- spear(s$score, h$meta$age)
  data.frame(cohort=g, n=ncol(h$expr), genes=s$n_up[1]+s$n_dn[1], rho=sp[1], p=sp[2]) }))
print(plat); print(sc_lvl); wr(as.data.frame(plat), "R1_platform_stratified.csv"); wr(sc_lvl, "R1_platform_array_signature_score.csv")
cat("RNA-seq-only meta: genes FDR<0.05 =", sum(rq$fdr<0.05), "\n")

## 5. Probe collapsing: mean-of-probes instead of highest-mean probe in array cohorts
alt <- harm; for (g in names(alt)) if (!is.null(alt[[g]]$expr_meanprobe)) alt[[g]]$expr <- alt[[g]]$expr_meanprobe
ea <- effects_all(alt); Ma <- zse(ea); r_alt <- meta_genes(Ma$Z, Ma$SE, 5, genes=intersect(rob$gene, rownames(Ma$Z)))
pa <- merge(rob[, c("gene","z","pval")], r_alt[, c("gene","z","pval","sign_consist")], by="gene", suffixes=c("_maxmean","_meanprobe"))
probe <- data.frame(genes=nrow(pa), same_direction=sum(sign(pa$z_maxmean)==sign(pa$z_meanprobe)), nominal_p05=sum(pa$pval_meanprobe<0.05),
                    r_effect=cor(pa$z_maxmean, pa$z_meanprobe), median_abs_r_ratio=median(abs(tanh(pa$z_meanprobe))/abs(tanh(pa$z_maxmean))))
print(probe); wr(probe, "R1_probe_collapsing.csv"); wr(pa, "R1_probe_collapsing_by_gene.csv")
saveRDS(list(r7=r7, rb7=rb7, ra=ra, sa=sa, rq=rq, th=th, dg=dg), file.path(D, "harmonized", paste0("sensitivity_", KSET, ".rds")))
