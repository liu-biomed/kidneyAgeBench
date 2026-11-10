# Exploratory pilots in the corrected discovery set to triage candidate new questions.
# Anything pursued must be locked here and then tested in independent data (TCGA / CPTAC / own tissue).
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
suppressMessages({library(msigdbr); library(mgcv); library(limma)})
b <- load_bundle(); harm <- b$harm; P <- readRDS(MP); res <- P$res
H <- as.data.frame(msigdbr(species="Homo sapiens", collection="H")); hl <- split(H$gene_symbol, H$gs_name)
sets <- list(
  TLS12_chemokine=c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13"),
  B_cell_TLS=c("MS4A1","CD79A","CD79B","CD19","CR2","FCRL5","CXCR5","LAMP3","LTB","CD22"),
  plasma_cell=c("JCHAIN","MZB1","IGHG1","IGKC","XBP1","TNFRSF17","SDC1","IGHA1"),
  T_cell=c("CD3D","CD3E","CD2","CD8A","GZMK","CD247","TRAC","IL7R"),
  macrophage=c("CD68","CD163","C1QA","C1QB","C1QC","MS4A7","CSF1R","FCGR3A"),
  interferon_gamma=hl$HALLMARK_INTERFERON_GAMMA_RESPONSE, oxphos=hl$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
  fatty_acid=hl$HALLMARK_FATTY_ACID_METABOLISM, PT_identity=c("LRP2","CUBN","SLC34A1","SLC5A12","SLC22A6","SLC22A8","HNF4A","ALDOB","PCK1","G6PC1"),
  fibrosis_ECM=c("COL1A1","COL1A2","COL3A1","COL6A3","FN1","POSTN","DCN","LUM","FBLN1","MGP"),
  injury_VCAM1_PT=c("VCAM1","HAVCR1","LCN2","SPP1","CD44","SOX9","KRT8","KRT18","CCL2","TNC"),
  SenMayo=unique(as.data.frame(msigdbr(species="Homo sapiens", collection="C2"))$gene_symbol[as.data.frame(msigdbr(species="Homo sapiens", collection="C2"))$gs_name=="SAUL_SEN_MAYO"]))
mz <- function(e, g){ g <- intersect(g, rownames(e)); z <- t(scale(t(e[g,,drop=FALSE]))); z[!is.finite(z)] <- NA; colMeans(z, na.rm=TRUE) }
dec <- read.csv(file.path(REV, "..", "tables", "deconv_proportions_bisque_class.csv"))
pc <- function(x, y, C) cor(residuals(lm(rank(x) ~ ., data=C)), residuals(lm(rank(y) ~ ., data=C)))
## 1. Programme-level age associations, raw and composition-adjusted
long <- list(); A <- list()
for (g in names(harm)) { h <- harm[[g]]; a <- h$meta$age; sx <- if (g %in% NO_SEX) NULL else factor(h$meta$sex)
  dd <- dec[match(first_gsm(h$meta$gsm), dec$gsm), c("epithelial","immune","stroma","endothelial")]
  for (s in names(sets)) { v <- mz(h$expr, sets[[s]])
    C0 <- if (is.null(sx)) data.frame(k=rep(0,length(a)))[, 0, drop=FALSE] else data.frame(sex=sx)
    r0 <- if (ncol(C0)) pc(a, v, C0) else cor(a, v, method="s")
    r1 <- if (all(!is.na(dd))) pc(a, v, cbind(C0, dd)) else NA
    A[[length(A)+1]] <- data.frame(cohort=g, n=length(a), set=s, genes=length(intersect(sets[[s]], rownames(h$expr))), rho=r0, rho_comp=r1)
    long[[length(long)+1]] <- data.frame(cohort=g, age=a, set=s, z=as.numeric(scale(v))) } }
A <- do.call(rbind, A)
prog <- do.call(rbind, lapply(names(sets), function(s){ d <- A[A$set==s,]; p0 <- pool_rho(d$rho, d$n); dc <- d[is.finite(d$rho_comp),]; p1 <- pool_rho(dc$rho_comp, dc$n-5)
  data.frame(set=s, genes_median=median(d$genes), rho=p0$rho, p=p0$p, I2=p0$I2, rho_comp=p1$rho, p_comp=p1$p, I2_comp=p1$I2) }))
prog <- prog[order(-abs(prog$rho_comp)),]; print(format(prog, digits=2), row.names=FALSE); wr(prog, "P1_programme_age_composition.csv")
## 2. Nonlinear trajectories: within-cohort z-scored programme ~ s(age) + cohort; compare with linear; locate acceleration
L <- do.call(rbind, long); L$cohort <- factor(L$cohort)
nl <- do.call(rbind, lapply(c("TLS12_chemokine","B_cell_TLS","plasma_cell","T_cell","macrophage","interferon_gamma","oxphos","PT_identity","fibrosis_ECM","SenMayo"), function(s){
  d <- L[L$set==s & L$age >= 20 & L$age <= 85,]; f1 <- lm(z ~ age + cohort, d); f2 <- gam(z ~ s(age, k=5) + cohort, data=d, method="REML")
  ag <- seq(25, 80, 1); nd <- data.frame(age=ag, cohort=levels(d$cohort)[1]); pr <- predict(f2, nd); sl <- diff(pr)
  # piecewise: best single breakpoint 35-75
  bp <- sapply(35:75, function(k) AIC(lm(z ~ age + pmax(age-k,0) + cohort, d))); k <- (35:75)[which.min(bp)]; fb <- lm(z ~ age + pmax(age-k,0) + cohort, d)
  data.frame(set=s, n=nrow(d), edf=summary(f2)$edf, p_nonlinear=anova(f1, gam(z ~ age + s(age, k=5, m=c(2,0)) + cohort, data=d, method="REML"), test="F")[2,"Pr(>F)"],
             AIC_linear=AIC(f1), AIC_gam=AIC(f2), breakpoint=k, slope_before=coef(fb)["age"]*10, slope_after=(coef(fb)["age"]+coef(fb)[3])*10, p_slope_change=summary(fb)$coef[3,4],
             max_slope_age=ag[which.max(abs(sl))]) }))
print(format(nl, digits=2), row.names=FALSE); wr(nl, "P2_nonlinear_trajectories.csv")
## 3. Ageing vs CKD: is the disease transcriptome an amplified ageing axis?
cl <- readRDS(file.path(D,"harmonized","clinical_bundle.rds"))
ck <- do.call(rbind, lapply(names(cl), function(co){ e <- cl[[co]]$expr; m <- cl[[co]]$meta; grp <- factor(m$grp, levels=c("Reference","Disease"))
  fit <- eBayes(lmFit(e, model.matrix(~grp))); tt <- topTable(fit, coef=2, number=Inf); tt$gene <- rownames(tt)
  mm <- merge(res[, c("gene","z","fdr")], tt[, c("gene","t")], by="gene"); sig <- mm[mm$fdr < 0.05,]
  data.frame(cohort=co, genes=nrow(mm), rho_genomewide=cor(mm$z, mm$t, method="s"), n_age_sig=nrow(sig),
             pct_age_sig_same_dir_disease=100*mean(sign(sig$z)==sign(sig$t)), rho_age_sig=cor(sig$z, sig$t, method="s"),
             up_shared=sum(sig$z>0 & sig$t>2), up_age_specific=sum(sig$z>0 & abs(sig$t)<1), dn_shared=sum(sig$z<0 & sig$t< -2), dn_age_specific=sum(sig$z<0 & abs(sig$t)<1)) }))
print(format(ck, digits=2), row.names=FALSE); wr(ck, "P3_ageing_vs_CKD_axis.csv")
