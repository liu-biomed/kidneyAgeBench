# Methodological and clinical upgrades for the Genome Medicine revision (manuscript v6):
#  1. Gene-level composition-adjusted benchmark (Z_comp) across all 48 sets and pooled gene sources (with and without CPTAC)
#  2. Inter-gene correlation (mean pairwise Spearman rho) and camera-style variance inflation factor (VIF = 1 + (k-1)*rho_bar)
#  3. Gene-level meta-analysis pooling sensitivity: Unweighted Stouffer vs Inverse-Variance Fixed-Effects (Liptak) vs DerSimonian-Laird Random-Effects
#  4. Reference-based single-nucleus deconvolution concordance (KPMP atlas: MuSiC, Bisque class/subclass, NNLS vs 5-compartment marker proxies in G6, n=414)
#  5. Direct clinical evaluation in three clinical settings (n = 428 clinical biopsies):
#     (a) ERCB CKD cohort (GSE104954, n = 190 across 9 CKD etiologies vs living donors)
#     (b) Tubulointerstitial RNA-seq CKD cohort with continuous eGFR (GSE175759, n = 87)
#     (c) Pre-implantation deceased-donor kidney biopsies with post-transplant graft function outcome (GSE147451, n = 151)

suppressMessages({library(data.table); library(metafor); library(ggplot2); library(patchwork)})
REV <- path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030")
source(file.path(REV, "scripts", "04_locked_definitions.R"))
source(file.path(REV, "scripts", "10_cptac_helpers.R"))
source(file.path(REV, "scripts", "40_fig_common.R"))
load(file.path(REV, "kidneyAgeBench/data/kab_library.rda"))
load(file.path(REV, "kidneyAgeBench/data/kab_bins.rda"))
load(file.path(REV, "kidneyAgeBench/data/kab_gene_z.rda"))
O <- file.path(REV, "benchmark", "results")

## ---- 1. Load normal kidney datasets and compute composition-adjusted gene-level age z ----
D <- list()
cb <- readRDS(path.expand("~/kidney_aging_data/harmonized/harmonized_bundle_corrected.rds"))
for (g in names(cb$harm)) {
  h <- cb$harm[[g]]
  D[[paste0("G6_", g)]] <- list(e=h$expr, age=h$meta$age, base=if (all(is.na(h$meta$sex))) NULL else data.frame(sex=factor(h$meta$sex)), grp="G6")
}
x <- fread(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_tpm.tsv"), data.table=FALSE)
tm <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_normal_meta.csv"))
sx <- read.csv(path.expand("~/kidney_aging_data/validation/tcga/tcga_case_sex.csv"))
tm$sex <- sx$sex[match(tm$patient, sx$patient)]
tm$age <- ifelse(is.finite(tm$age_days), tm$age_days/365.25, tm$age_at_index)
tm <- tm[order(tm$sample), ]; tm <- tm[!duplicated(tm$patient) & is.finite(tm$age) & !is.na(tm$sex), ]
e <- log2(as.matrix(x[, tm$file_id]) + 1); o <- order(rowMeans(e), decreasing=TRUE); e <- e[o,]; gg <- x$gene_name[o]
e <- e[!duplicated(gg),]; rownames(e) <- gg[!duplicated(gg)]; e <- e[rowMeans(e > 0) >= 0.2, ]
D$TCGA <- list(e=e, age=tm$age, base=data.frame(sex=factor(tm$sex), project=factor(tm$project)), grp="TCGA")

V <- path.expand("~/kidney_aging_data/validation/cptac")
mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
tr <- cptac_rd(V, "bcm-ccrcc-transcriptomics"); tr <- cbind(gene=tr[[1]], tr[, grep("_A$", names(tr))])
R <- cptac_sym(tr, e2s); colnames(R) <- sub("_A$", "", colnames(R)); R <- R[rowMeans(R > 0) >= 0.2, ]
cs <- intersect(colnames(R), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")]); mm <- C[match(cs, C$case), ]
D$CPTAC <- list(e=R[, cs], age=mm$age, base=data.frame(sex=factor(mm$sex)), grp="CPTAC")
rm(x, e, R, tr); invisible(gc())

gene_z_cov <- function(d, add_comp=FALSE) {
  px <- if (add_comp) as.data.frame(lapply(LOCK$proxy, function(g) mz(d$e, g))) else NULL
  covs <- if (is.null(d$base) && is.null(px)) NULL else if (is.null(d$base)) px else if (is.null(px)) d$base else cbind(d$base, px)
  if (!is.null(covs)) covs <- covs[, sapply(covs, function(v) all(is.finite(v) | is.factor(v))), drop=FALSE]
  X <- if (is.null(covs) || ncol(covs)==0) matrix(1, length(d$age)) else model.matrix(~., data=covs)
  q <- ncol(X) - 1
  Hm <- X %*% solve(crossprod(X)) %*% t(X); a <- rank(d$age, na.last="keep"); a <- a - Hm %*% a
  Rk <- t(apply(d$e, 1, rank, na.last="keep")); Rk <- Rk - t(Hm %*% t(Rk)); r <- as.numeric((Rk/sqrt(rowSums(Rk^2))) %*% (a/sqrt(sum(a^2))))
  setNames(atanh(pmin(pmax(r, -.9999), .9999)) * sqrt(length(d$age) - q - 3), rownames(d$e))
}

Z_base <- lapply(D, gene_z_cov, add_comp=FALSE)
Z_comp <- lapply(D, gene_z_cov, add_comp=TRUE)
grp <- sapply(D, `[[`, "grp")
min_k <- function(n) max(2, floor(n / 2))
pooled <- function(Zlist, grps) {
  use <- names(Zlist)[grp %in% grps]
  g <- Reduce(union, lapply(Zlist[use], names)); M <- sapply(Zlist[use], function(z) z[g]); rownames(M) <- g
  k <- rowSums(is.finite(M)); s <- rowSums(M, na.rm=TRUE)/sqrt(pmax(k,1)); s[k < min_k(length(use))] <- NA; s
}
eval_bins <- function(use) {
  P <- lapply(D[use], function(d) { m <- d$e; mu <- rowMeans(m); sdv <- apply(m, 1, sd); cbind(mu = rank(mu)/length(mu), sd = rank(sdv)/length(sdv)) })
  g <- Reduce(union, lapply(P, rownames)); mu <- sapply(P, function(p) p[, "mu"][g]); sdv <- sapply(P, function(p) p[, "sd"][g]); rownames(mu) <- rownames(sdv) <- g
  keep <- rowSums(is.finite(mu)) >= min_k(length(use)); m1 <- rowMeans(mu[keep, ], na.rm=TRUE); s1 <- rowMeans(sdv[keep, ], na.rm=TRUE)
  setNames(paste(cut(rank(m1), 10, labels=FALSE), cut(rank(s1), 10, labels=FALSE)), names(m1))
}
bins_cop <- eval_bins(names(D)[grp %in% c("G6", "TCGA")])
leak <- function(id) if (grepl("^A1[12]_", id)) "G6" else if (grepl("^A0[56]_Rowland", id)) "TCGA" else character()

## ---- 2. Inter-gene correlation (camera VIF) and composition-adjusted 48-set benchmark ----
sigs <- split(kab_library, kab_library$sig_id)
igc_map <- sapply(sigs, function(s) {
  g <- unique(s$gene)
  rb <- sapply(D, function(d) {
    gg <- intersect(g, rownames(d$e)); if (length(gg) < 2) return(NA_real_)
    Cmat <- cor(t(d$e[gg, , drop=FALSE]), method="s")
    mean(Cmat[upper.tri(Cmat)], na.rm=TRUE)
  })
  mean(rb, na.rm=TRUE)
})

set.seed(20261101)
run_48_comp <- function(Zlist, grps, bins, cfg_name) {
  res <- do.call(rbind, lapply(sigs, function(s) {
    id <- s$sig_id[1]; use_g <- setdiff(grps, leak(id))
    if (length(names(Zlist)[grp %in% use_g]) < 2) return(NULL)
    pz <- pooled(Zlist, use_g)
    ok <- intersect(names(pz)[is.finite(pz)], names(bins)); g <- intersect(unique(s$gene), ok)
    if (length(g) < 2) return(NULL)
    pool <- split(ok, bins[ok]); obs <- mean(abs(pz[g]))
    idx <- vapply(g, function(x) { p <- pool[[bins[[x]]]]; p[sample.int(length(p), 2000, TRUE)] }, character(2000))
    Zn <- matrix(pz[idx], nrow=2000)
    nul <- rowMeans(abs(Zn))
    rb <- igc_map[[id]]; vif <- 1 + (length(g) - 1) * max(0, rb)
    z_emp <- (obs - mean(nul)) / sd(nul)
    z_vif <- z_emp / sqrt(vif)
    p_vif <- pnorm(z_vif, lower.tail=FALSE)
    # Matched-null directional test for signed sets
    dir_v <- s$dir[match(g, s$gene)]; sg_mask <- !is.na(dir_v) & dir_v != 0
    if (any(sg_mask)) {
      D_obs <- mean(dir_v[sg_mask] * pz[g[sg_mask]])
      D_nul <- as.numeric(Zn[, sg_mask, drop=FALSE] %*% dir_v[sg_mask]) / sum(sg_mask)
      D_z <- (D_obs - mean(D_nul)) / sd(D_nul)
      p_dir_matched <- (1 + sum(D_nul >= D_obs)) / 2001
    } else { D_z <- NA_real_; p_dir_matched <- NA_real_ }
    data.frame(config=cfg_name, sig_id=id, class=s$class[1], n=length(g),
               ER=obs/mean(nul), emp_p=(1+sum(nul >= obs))/2001,
               mean_intergene_rho=rb, camera_VIF=vif, p_vif=p_vif,
               D_excess_z=D_z, p_dir_matched=p_dir_matched)
  }))
  res$fdr <- p.adjust(res$emp_p, "BH")
  res$fdr_vif <- p.adjust(res$p_vif, "BH")
  res
}

comp_48 <- rbind(
  run_48_comp(Z_base, c("G6","TCGA","CPTAC"), kab_bins, "base_locked"),
  run_48_comp(Z_comp, c("G6","TCGA","CPTAC"), kab_bins, "comp_locked"),
  run_48_comp(Z_base, c("G6","TCGA"), bins_cop, "base_coprimary"),
  run_48_comp(Z_comp, c("G6","TCGA"), bins_cop, "comp_coprimary")
)
write.csv(comp_48, file.path(O, "REVISION_gene_level_composition_adjusted.csv"), row.names=FALSE)

## ---- 3. Pooled gene sources in Unadjusted Age, Composition-Adjusted Age, and CKD Biopsies ----
cl <- readRDS(path.expand("~/kidney_aging_data/harmonized/clinical_bundle.rds"))
z_vec <- function(e, y, C=NULL) {
  ok <- is.finite(y); e <- e[, ok]; y <- y[ok]
  X <- if (is.null(C)) matrix(1, length(y)) else model.matrix(~., data=C[ok, , drop=FALSE]); q <- ncol(X) - 1
  Hm <- X %*% solve(crossprod(X)) %*% t(X); a <- rank(y, na.last="keep"); a <- a - Hm %*% a
  Rk <- t(apply(e, 1, rank, na.last="keep")); Rk <- Rk - t(Hm %*% t(Rk)); r <- as.numeric((Rk/sqrt(rowSums(Rk^2))) %*% (a/sqrt(sum(a^2))))
  setNames(atanh(pmin(pmax(r, -.9999), .9999)) * sqrt(length(y) - q - 3), rownames(e))
}
keep104_pool <- cl$GSE104954$meta$diagnosis != "Tumor nephrectomy"
e104_pool <- cl$GSE104954$expr[, keep104_pool]; m104_pool <- cl$GSE104954$meta[keep104_pool, ]
z_104 <- z_vec(e104_pool, as.numeric(m104_pool$grp == "Disease"), data.frame(plat=factor(m104_pool$platform)))
z_175_ckd <- z_vec(cl$GSE175759$expr, as.numeric(cl$GSE175759$meta$grp == "Disease"))
g_ckd <- intersect(names(z_104), names(z_175_ckd))
z_ckd_pool <- (z_104[g_ckd] + z_175_ckd[g_ckd]) / sqrt(2)
mu104 <- rowMeans(e104_pool[g_ckd, ]); sd104 <- apply(e104_pool[g_ckd, ], 1, sd)
mu175 <- rowMeans(cl$GSE175759$expr[g_ckd, ]); sd175 <- apply(cl$GSE175759$expr[g_ckd, ], 1, sd)
bins_ckd <- setNames(paste(cut(rank((rank(mu104)+rank(mu175))/2), 10, labels=FALSE), cut(rank((rank(sd104)+rank(sd175))/2), 10, labels=FALSE)), g_ckd)

lib <- kab_library
raag_all <- unique(lib$gene[lib$sig_id == "A08_Perco2019_RAAG634"])
raag_ckd <- unique(lib$gene[lib$sig_id == "A09_Perco2019_RAAG_CKD"])
pools <- list(
  "Class A (external age sets)" = unique(lib$gene[lib$class == "A" & !grepl("^A1[12]_|^A0[56]_", lib$sig_id)]),
  "RAAG, CKD-concordant" = raag_ckd,
  "RAAG, non-CKD" = setdiff(raag_all, raag_ckd),
  "Class C ('ageing' disease sets)" = unique(lib$gene[lib$class == "C"]),
  "CellAge" = unique(lib$gene[lib$sig_id == "B04_CellAge"]),
  "GenAge" = unique(lib$gene[lib$sig_id == "B03_GenAge"])
)

eval_pool_list <- list(
  "Age (unadjusted, locked)" = list(z = pooled(Z_base, c("G6","TCGA","CPTAC")), b = kab_bins),
  "Age (composition-adjusted, locked)" = list(z = pooled(Z_comp, c("G6","TCGA","CPTAC")), b = kab_bins),
  "Age (unadjusted, co-primary)" = list(z = pooled(Z_base, c("G6","TCGA")), b = bins_cop),
  "Age (composition-adjusted, co-primary)" = list(z = pooled(Z_comp, c("G6","TCGA")), b = bins_cop),
  "CKD vs reference (GSE104954 + GSE175759)" = list(z = z_ckd_pool, b = bins_ckd)
)

set.seed(20261101)
pool_df <- do.call(rbind, lapply(names(eval_pool_list), function(an) {
  zv <- eval_pool_list[[an]]$z; bv <- eval_pool_list[[an]]$b
  ok <- intersect(names(zv)[is.finite(zv)], names(bv)); pl <- split(ok, bv[ok])
  do.call(rbind, lapply(names(pools), function(pn) {
    g <- intersect(pools[[pn]], ok); obs <- mean(abs(zv[g]))
    nul <- rowMeans(vapply(g, function(x) { p <- pl[[bv[[x]]]]; abs(zv[p[sample.int(length(p), 2000, TRUE)]]) }, numeric(2000)))
    data.frame(analysis=an, pool=pn, n=length(g), ER=obs/mean(nul),
               lo=obs/quantile(nul, .975), hi=obs/quantile(nul, .025), emp_p=(1+sum(nul >= obs))/2001)
  }))
}))
write.csv(pool_df, file.path(O, "REVISION_gene_pools_comp_and_ckd.csv"), row.names=FALSE)

## ---- 4. Meta-analysis pooling sensitivity (Stouffer vs Fixed-Effects vs Random-Effects DL, with leakage control) ----
dfs <- sapply(D, function(d) {
  q <- if (is.null(d$base)) 0 else ncol(model.matrix(~., data=d$base)) - 1
  length(d$age) - q - 3
})
pool_three_methods <- function(excl) {
  use_m <- names(Z_base)[grp %in% setdiff(c("G6","TCGA","CPTAC"), excl)]
  g_all <- Reduce(union, lapply(Z_base[use_m], names))
  M <- sapply(Z_base[use_m], function(z) z[g_all]); rownames(M) <- g_all
  keep_g <- rowSums(is.finite(M)) >= min_k(length(use_m))
  M <- M[keep_g, , drop=FALSE]
  z_st <- rowSums(M, na.rm=TRUE) / sqrt(rowSums(is.finite(M)))
  W_mat <- matrix(sqrt(dfs[use_m]), nrow=nrow(M), ncol=ncol(M), byrow=TRUE); W_mat[!is.finite(M)] <- NA
  z_fe <- rowSums(M * W_mat, na.rm=TRUE) / sqrt(rowSums(W_mat^2, na.rm=TRUE))
  yi <- M / W_mat; vi <- 1 / (W_mat^2); wi <- 1 / vi
  y_fe <- rowSums(wi * yi, na.rm=TRUE) / rowSums(wi, na.rm=TRUE)
  Q <- rowSums(wi * (yi - y_fe)^2, na.rm=TRUE)
  c_dl <- rowSums(wi, na.rm=TRUE) - rowSums(wi^2, na.rm=TRUE) / rowSums(wi, na.rm=TRUE)
  tau2 <- pmax(0, (Q - (rowSums(is.finite(yi)) - 1)) / c_dl)
  wi_re <- 1 / (vi + tau2)
  z_re <- (rowSums(wi_re * yi, na.rm=TRUE) / rowSums(wi_re, na.rm=TRUE)) / sqrt(1 / rowSums(wi_re, na.rm=TRUE))
  list(st=z_st, fe=z_fe, re=z_re)
}
PM_list <- list(none=pool_three_methods(character()), G6=pool_three_methods("G6"), TCGA=pool_three_methods("TCGA"))

set.seed(20261101)
meta_sens <- do.call(rbind, lapply(sigs, function(s) {
  id <- s$sig_id[1]; lk <- leak(id)
  pm <- PM_list[[if (length(lk)) lk else "none"]]
  ok_m <- intersect(names(pm$st)[is.finite(pm$st)], names(kab_bins))
  pool_m <- split(ok_m, kab_bins[ok_m])
  g <- intersect(unique(s$gene), ok_m); if (length(g) < 2) return(NULL)
  idx <- vapply(g, function(x) { p <- pool_m[[kab_bins[[x]]]]; p[sample.int(length(p), 2000, TRUE)] }, character(2000))
  f_er <- function(zv) { nul <- rowMeans(matrix(abs(zv[idx]), nrow=2000)); obs <- mean(abs(zv[g])); c(obs/mean(nul), (1+sum(nul >= obs))/2001) }
  r1 <- f_er(pm$st); r2 <- f_er(pm$fe); r3 <- f_er(pm$re)
  data.frame(sig_id=id, class=s$class[1], n=length(g),
             ER_Stouffer=r1[1], p_Stouffer=r1[2], ER_FE_weighted=r2[1], p_FE_weighted=r2[2], ER_RE_DL=r3[1], p_RE_DL=r3[2])
}))
write.csv(meta_sens, file.path(O, "REVISION_meta_pooling_sensitivity.csv"), row.names=FALSE)

## ---- 5. KPMP single-nucleus reference deconvolution vs 5-marker composition adjustment in G6 ----
bisq <- read.csv(file.path(REV, "..", "tables", "deconv_proportions_bisque_class.csv"))
bisq_sub <- read.csv(file.path(REV, "..", "tables", "deconv_proportions_bisque_subclass.csv"))
musc <- read.csv(file.path(REV, "..", "tables", "deconv_proportions_music_class.csv"))
nnls <- read.csv(file.path(REV, "..", "tables", "deconv_proportions_nnls_class.csv"))
first_gsm <- function(x) sub(";.*", "", x)
signed_ids <- names(sigs)[sapply(sigs, function(s) any(s$dir != 0))]
score_sig <- function(e, s) {
  up <- s$gene[s$dir > 0]; dn <- s$gene[s$dir < 0]
  u <- mz(e, up, 2); d <- mz(e, dn, 2)
  if (all(is.na(u)) && all(is.na(d))) return(NULL)
  if (all(is.na(d))) u else if (all(is.na(u))) -d else u - d
}
pool_g6 <- function(r, n, q) {
  ok <- is.finite(r)
  f <- tryCatch(rma(yi=atanh(r[ok]), sei=1/sqrt(n[ok]-q[ok]-3), method="REML"),
                error=function(e) rma(yi=atanh(r[ok]), sei=1/sqrt(n[ok]-q[ok]-3), method="DL"))
  tanh(as.numeric(f$b))
}
use_c <- intersect(names(cb$harm), unique(bisq$cohort))
deconv_cmp <- do.call(rbind, lapply(setdiff(signed_ids, c("A11_KAS94_original", "A12_KAS82_corrected")), function(id) {
  s <- sigs[[id]]
  out_c <- do.call(rbind, lapply(use_c, function(g) {
    h <- cb$harm[[g]]; v <- score_sig(h$expr, s); if (is.null(v)) return(NULL)
    base <- if (all(is.na(h$meta$sex))) NULL else data.frame(sex=factor(h$meta$sex))
    px <- as.data.frame(lapply(LOCK$proxy, function(gg) mz(h$expr, gg)))
    fg <- first_gsm(h$meta$gsm)
    b_cov <- bisq[match(fg, bisq$gsm), c("immune","stroma","endothelial"), drop=FALSE]
    bs_cov <- bisq_sub[match(fg, bisq_sub$gsm), c("PT","DCT","IMM","FIB","EC"), drop=FALSE]
    m_cov <- musc[match(fg, musc$gsm), c("immune","stroma","endothelial"), drop=FALSE]
    n_cov <- nnls[match(fg, nnls$gsm), c("immune","stroma","endothelial"), drop=FALSE]
    c_add <- function(b, x) if (is.null(b)) as.data.frame(x) else cbind(b, as.data.frame(x))
    data.frame(
      n = length(h$meta$age),
      q_base = if (is.null(base)) 0 else 1,
      q_px = (if (is.null(base)) 0 else 1) + 5,
      q_dec = (if (is.null(base)) 0 else 1) + 3,
      r_base = pspear(h$meta$age, v, base)$rho,
      r_px = pspear(h$meta$age, v, c_add(base, px))$rho,
      r_bisq_sub = pspear(h$meta$age, v, c_add(base, bs_cov))$rho,
      r_bisq = pspear(h$meta$age, v, c_add(base, b_cov))$rho,
      r_musc = pspear(h$meta$age, v, c_add(base, m_cov))$rho,
      r_nnls = pspear(h$meta$age, v, c_add(base, n_cov))$rho
    )
  }))
  data.frame(
    sig_id = id,
    G6_unadjusted = pool_g6(out_c$r_base, out_c$n, out_c$q_base),
    G6_marker5 = pool_g6(out_c$r_px, out_c$n, out_c$q_px),
    G6_KPMP_MuSiC = pool_g6(out_c$r_musc, out_c$n, out_c$q_dec),
    G6_KPMP_Bisque_5sub = pool_g6(out_c$r_bisq_sub, out_c$n, out_c$q_px),
    G6_KPMP_Bisque_class = pool_g6(out_c$r_bisq, out_c$n, out_c$q_dec),
    G6_KPMP_NNLS = pool_g6(out_c$r_nnls, out_c$n, out_c$q_dec)
  )
}))
write.csv(deconv_cmp, file.path(O, "REVISION_deconv_method_concordance.csv"), row.names=FALSE)

## ---- 6. Direct clinical evaluation in CKD biopsies (GSE104954, GSE175759) and donor biopsies (GSE147451) ----
inj_genes <- LOCK$outcome$PT_injury
e104 <- cl$GSE104954$expr; m104 <- cl$GSE104954$meta
keep104 <- m104$diagnosis != "Tumor nephrectomy"
e104 <- e104[, keep104]; m104 <- m104[keep104, ]
inj104 <- mz(e104, inj_genes); px104 <- as.data.frame(lapply(LOCK$proxy, function(g) mz(e104, g)))
y104 <- as.numeric(m104$grp == "Disease")

e175 <- cl$GSE175759$expr; m175 <- cl$GSE175759$meta
inj175 <- mz(e175, inj_genes); px175 <- as.data.frame(lapply(LOCK$proxy, function(g) mz(e175, g)))
y175 <- as.numeric(m175$grp == "Disease")

h147 <- cb$harm$GSE147451; e147 <- h147$expr; m147 <- h147$meta
px147 <- as.data.frame(lapply(LOCK$proxy, function(g) mz(e147, g)))
ok147 <- !is.na(m147$gfr_group); y147 <- as.numeric(m147$gfr_group[ok147] == "poor")

clin_tab <- do.call(rbind, lapply(signed_ids, function(id) {
  # GSE104954
  v1 <- score_sig(e104, sigs[[id]]); v1z <- (v1 - mean(v1[y104==0])) / sd(v1[y104==0])
  f1_raw <- summary(lm(v1z ~ y104 + factor(m104$platform)))$coef["y104", ]
  f1_adj <- summary(lm(v1z ~ y104 + factor(m104$platform) + inj104 + ., data=px104))$coef["y104", ]
  r_inj104 <- cor(v1, inj104, method="s")
  # GSE175759
  v2 <- score_sig(e175, sigs[[id]]); v2z <- (v2 - mean(v2[y175==0])) / sd(v2[y175==0])
  f2_raw <- summary(lm(v2z ~ y175))$coef["y175", ]
  f2_adj <- summary(lm(v2z ~ y175 + inj175 + ., data=px175))$coef["y175", ]
  eg_raw <- pspear(m175$egfr, v2)
  eg_adj <- pspear(m175$egfr, v2, cbind(px175, inj=inj175))
  # GSE147451
  v3 <- score_sig(e147, sigs[[id]])
  tx_raw <- pspear(y147, v3[ok147])
  tx_age <- pspear(y147, v3[ok147], data.frame(age=m147$age[ok147]))
  tx_comp <- pspear(y147, v3[ok147], cbind(age=m147$age[ok147], px147[ok147, ]))
  data.frame(
    sig_id = id,
    ERCB_beta_raw = f1_raw[1], ERCB_se_raw = f1_raw[2], ERCB_p_raw = f1_raw[4],
    ERCB_beta_adj = f1_adj[1], ERCB_se_adj = f1_adj[2], ERCB_p_adj = f1_adj[4],
    ERCB_cor_PT_injury = r_inj104,
    GSE175_beta_raw = f2_raw[1], GSE175_p_raw = f2_raw[4],
    GSE175_beta_adj = f2_adj[1], GSE175_p_adj = f2_adj[4],
    GSE175_eGFR_rho_raw = eg_raw$rho, GSE175_eGFR_p_raw = eg_raw$p,
    GSE175_eGFR_rho_adj = eg_adj$rho, GSE175_eGFR_p_adj = eg_adj$p,
    TX_poorGFR_rho_raw = tx_raw$rho, TX_poorGFR_p_raw = tx_raw$p,
    TX_poorGFR_rho_adj_age = tx_age$rho, TX_poorGFR_p_adj_age = tx_age$p,
    TX_poorGFR_rho_adj_age_comp = tx_comp$rho, TX_poorGFR_p_adj_age_comp = tx_comp$p
  )
}))
write.csv(clin_tab, file.path(O, "REVISION_clinical_CKD_and_transplant.csv"), row.names=FALSE)

## ---- 7. Build Figure 8 (Clinical translation & gene-level composition resolution) ----
# Panel a: Pooled gene sources across Unadjusted Age, Composition-Adjusted Age, and CKD vs Reference
pa_df <- pool_df[pool_df$analysis %in% c("Age (unadjusted, locked)", "Age (composition-adjusted, locked)", "CKD vs reference (GSE104954 + GSE175759)"), ]
pa_df$an_lab <- factor(c(
  "Age (unadjusted, locked)" = "Normal kidney age (unadjusted)",
  "Age (composition-adjusted, locked)" = "Normal kidney age (composition-adjusted)",
  "CKD vs reference (GSE104954 + GSE175759)" = "CKD vs reference biopsies (n = 277)"
)[pa_df$analysis], levels=c("Normal kidney age (unadjusted)", "Normal kidney age (composition-adjusted)", "CKD vs reference biopsies (n = 277)"))
pa_df$pool <- factor(pa_df$pool, levels=rev(c("Class A (external age sets)", "RAAG, CKD-concordant", "RAAG, non-CKD", "Class C ('ageing' disease sets)", "CellAge", "GenAge")))

p8a <- ggplot(pa_df, aes(ER, pool, colour=an_lab, shape=an_lab)) +
  geom_vline(xintercept=1, linetype=2, linewidth=.3, colour="grey40") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=.4, orientation="y", position=position_dodge(width=.62)) +
  geom_point(size=1.9, stroke=.45, fill="white", position=position_dodge(width=.62)) +
  scale_colour_manual(values=c("Normal kidney age (unadjusted)"="#1D63B8", "Normal kidney age (composition-adjusted)"="#008B7D", "CKD vs reference biopsies (n = 277)"="#C87028"), name=NULL) +
  scale_shape_manual(values=c(16, 21, 17), name=NULL) +
  labs(x="Enrichment vs matched genes (95% null interval)", y=NULL) +
  theme(legend.position="top", legend.direction="vertical", legend.margin=margin(0,0,0,0))

# Panel b: ERCB 9 CKD etiologies vs Living donor ordered by PT injury
diag_map <- c(
  "Living donor" = "Living donor",
  "Thin membrande disease" = "TMD",
  "Minimal change disease" = "MCD",
  "Membranous glomerulonephropathy" = "MN",
  "IgA nephropathy" = "IgAN",
  "Hypertensive nephropathy" = "HTN",
  "Systemic lupus erythematosus" = "LN",
  "Focal segmental glomerulosclerosis" = "FSGS",
  "Diabetic nephropathy" = "DKD",
  "ANCA-associated vasculitis" = "ANCA"
)
keep_d <- m104$diagnosis %in% names(diag_map)
show_sigs <- c("A01_Rodwell2004"="Rodwell 2004", "A04_Melk2005_PAM"="Melk 2005 (PAM)", "A06_Rowland2019_rep19"="Rowland 2019 (rep)", "A12_KAS82_corrected"="KAS-82", "B01_deMagalhaes2009"="de Magalhães 2009", "B05_SenMayo"="SenMayo")
sc_m <- sapply(names(show_sigs), function(id) {
  v <- score_sig(e104, sigs[[id]]); (v - mean(v[m104$diagnosis=="Living donor"])) / sd(v[m104$diagnosis=="Living donor"])
})
pt_z <- (inj104 - mean(inj104[m104$diagnosis=="Living donor"])) / sd(inj104[m104$diagnosis=="Living donor"])
pb_raw <- cbind(data.frame(diag=diag_map[m104$diagnosis[keep_d]], PT_injury=pt_z[keep_d]), as.data.frame(sc_m[keep_d, ]))
pb_agg <- aggregate(. ~ diag, pb_raw, median)
ord_diag <- pb_agg$diag[order(pb_agg$PT_injury)]
pb_long <- reshape(pb_agg, varying=names(show_sigs), v.names="score_SD", timevar="sig_id", times=names(show_sigs), direction="long")
pb_long$diag <- factor(pb_long$diag, levels=ord_diag)
pb_long$sig_lab <- factor(show_sigs[pb_long$sig_id], levels=unname(show_sigs))
pb_long$cls <- ifelse(grepl("^A", pb_long$sig_id), "A  Kidney age-derived", "B  General ageing / senescence")

# End-of-line labels at ANCA (x = max PT_injury) and disease annotations along top
end_df <- do.call(rbind, lapply(split(pb_long, pb_long$sig_lab), function(d) {
  f <- lm(score_SD ~ PT_injury, data=d)
  mx <- max(pb_agg$PT_injury)
  data.frame(sig_lab=d$sig_lab[1], cls=d$cls[1], PT_injury=mx, score_SD=as.numeric(predict(f, newdata=data.frame(PT_injury=mx))))
}))
end_df <- end_df[order(end_df$score_SD), ]
min_gap <- 0.36
for (i in 2:nrow(end_df)) {
  if (end_df$score_SD[i] - end_df$score_SD[i - 1] < min_gap) {
    end_df$score_SD[i] <- end_df$score_SD[i - 1] + min_gap
  }
}
# Shift top cluster down slightly if needed so all stay neatly centered around their endpoints
end_df$y_lab <- end_df$score_SD - c(0, 0, 0.14, 0.14, 0.14, 0.14)
end_df$score_end <- sapply(split(pb_long, pb_long$sig_lab)[as.character(end_df$sig_lab)], function(d) {
  as.numeric(predict(lm(score_SD ~ PT_injury, data=d), newdata=data.frame(PT_injury=max(pb_agg$PT_injury))))
})
# Merge close x values for disease guide labels
diag_ticks <- data.frame(
  PT_injury = c(-0.28, 0.33, 0.68, 1.58, 1.89, 2.20, 2.82, 4.01),
  lab = c("Living\ndonor", "TMD /\nMCD", "MN", "IgAN /\nHTN", "LN", "FSGS", "DKD", "ANCA")
)

p8b <- ggplot(pb_long, aes(PT_injury, score_SD, group=sig_lab, colour=cls)) +
  geom_vline(data=diag_ticks, aes(xintercept=PT_injury), linetype=3, linewidth=.22, colour="grey78") +
  geom_hline(yintercept=0, linetype=2, linewidth=.25, colour="grey55") +
  geom_smooth(method="lm", se=FALSE, linewidth=.4) +
  geom_point(size=1.35, alpha=.85) +
  geom_text(data=diag_ticks, aes(x=PT_injury, y=5.25, label=lab), inherit.aes=FALSE, size=TX*0.85, vjust=1, lineheight=0.85, colour="grey30") +
  geom_segment(data=end_df, aes(x=PT_injury + 0.03, xend=PT_injury + 0.10, y=score_end, yend=y_lab), inherit.aes=FALSE, linewidth=.22, colour="grey60") +
  geom_text(data=end_df, aes(x=PT_injury + 0.12, y=y_lab, label=sig_lab), inherit.aes=FALSE, colour="grey15",
            hjust=0, size=TX*0.92, show.legend=FALSE) +
  scale_colour_manual(values=c("A  Kidney age-derived"="#1D63B8", "B  General ageing / senescence"="#C87028"), name=NULL) +
  scale_x_continuous(limits=c(-0.45, 5.05), breaks=0:4) +
  scale_y_continuous(limits=c(-0.8, 5.40), breaks=0:4) +
  labs(x="Median proximal tubular injury score by CKD etiology in ERCB (SD vs living donor, n = 190)",
       y="Median signature score\n(SD vs living donor)") +
  theme(legend.position="top", legend.direction="horizontal", legend.margin=margin(0,0,0,0))

# Panel c: Attenuation of apparent CKD "accelerated ageing" in ERCB (GSE104954) after adjusting for composition + PT injury
pc_df <- clin_tab
pc_df$lab <- factor(F$lab[match(pc_df$sig_id, F$sig_id)], levels=F$lab[match(pc_df$sig_id[order(pc_df$ERCB_beta_raw)], F$sig_id)])
pc_long <- rbind(
  data.frame(lab=pc_df$lab, model="Unadjusted", beta=pc_df$ERCB_beta_raw, se=pc_df$ERCB_se_raw),
  data.frame(lab=pc_df$lab, model="Adjusted for composition + PT injury", beta=pc_df$ERCB_beta_adj, se=pc_df$ERCB_se_adj)
)
pc_long$model <- factor(pc_long$model, levels=c("Unadjusted", "Adjusted for composition + PT injury"))

p8c <- ggplot(pc_long, aes(beta, lab, shape=model, fill=model)) +
  geom_vline(xintercept=0, linetype=2, linewidth=.3, colour="grey40") +
  geom_line(aes(group=lab), linewidth=.35, colour="grey55") +
  geom_point(size=1.8, stroke=.45, colour="#1F5A8A") +
  scale_shape_manual(values=c(16, 21), name=NULL) +
  scale_fill_manual(values=c("#1F5A8A", "white"), name=NULL) +
  labs(x="Apparent ageing shift in CKD vs living donor (SD, GSE104954)", y=NULL) +
  theme(legend.position="top", legend.direction="vertical", legend.margin=margin(0,0,0,0))

# Panel d: Pre-implantation donor biopsies (GSE147451, n=151): predicting poor post-transplant GFR
pd_df <- clin_tab
pd_df$lab <- factor(F$lab[match(pd_df$sig_id, F$sig_id)], levels=F$lab[match(pd_df$sig_id[order(pd_df$TX_poorGFR_rho_raw)], F$sig_id)])
pd_long <- rbind(
  data.frame(lab=pd_df$lab, model="Unadjusted", rho=pd_df$TX_poorGFR_rho_raw, sig=pd_df$TX_poorGFR_p_raw < 0.05),
  data.frame(lab=pd_df$lab, model="Adjusted for donor age", rho=pd_df$TX_poorGFR_rho_adj_age, sig=pd_df$TX_poorGFR_p_adj_age < 0.05),
  data.frame(lab=pd_df$lab, model="Adjusted for donor age + composition", rho=pd_df$TX_poorGFR_rho_adj_age_comp, sig=pd_df$TX_poorGFR_p_adj_age_comp < 0.05)
)
pd_long$model <- factor(pd_long$model, levels=c("Unadjusted", "Adjusted for donor age", "Adjusted for donor age + composition"))

p8d <- ggplot(pd_long, aes(rho, lab, colour=model, shape=model)) +
  geom_vline(xintercept=0, linetype=2, linewidth=.3, colour="grey40") +
  geom_vline(xintercept=0.301, linetype=3, linewidth=.35, colour="#B2182B") +
  geom_line(aes(group=lab), linewidth=.3, colour="grey70") +
  geom_point(size=1.7, stroke=.45, fill="white") +
  annotate("text", x=0.295, y=2.2, label="Donor age\nalone (ρ = 0.30)", hjust=1, size=TX*0.9, colour="grey20") +
  scale_colour_manual(values=c("#1D63B8", "#C87028", "#008B7D"), name=NULL) +
  scale_shape_manual(values=c(16, 17, 21), name=NULL) +
  labs(x="Correlation with poor graft function (ρ, GSE147451)", y=NULL) +
  theme(legend.position="top", legend.direction="vertical", legend.margin=margin(0,0,0,0))

# Panel e: 5-marker composition adjustment vs KPMP single-nucleus reference deconvolution (MuSiC and Bisque) in G6 (n=414)
pe_df <- rbind(
  data.frame(sig_id=deconv_cmp$sig_id, marker5=deconv_cmp$G6_marker5, deconv=deconv_cmp$G6_KPMP_MuSiC, method="MuSiC (KPMP snRNA-seq, ρ = 0.95)"),
  data.frame(sig_id=deconv_cmp$sig_id, marker5=deconv_cmp$G6_marker5, deconv=deconv_cmp$G6_KPMP_Bisque_5sub, method="Bisque 5-subclass (KPMP snRNA-seq, ρ = 0.87)")
)
p8e <- ggplot(pe_df, aes(marker5, deconv, colour=method, shape=method)) +
  geom_abline(slope=1, intercept=0, linetype=2, linewidth=.3, colour="grey50") +
  geom_hline(yintercept=0, linetype=3, linewidth=.25, colour="grey65") +
  geom_vline(xintercept=0, linetype=3, linewidth=.25, colour="grey65") +
  geom_point(size=1.9, stroke=.45, fill="white") +
  scale_colour_manual(values=c("#1D63B8", "#C87028"), name=NULL) +
  scale_shape_manual(values=c(16, 21), name=NULL) +
  labs(x="G6 age ρ (5-compartment marker adjusted)",
       y="G6 age ρ (KPMP snRNA-seq\ndeconvolution adjusted)") +
  theme(legend.position="top", legend.direction="vertical", legend.margin=margin(0,0,0,0))

fig8 <- (p8a + p8e + plot_layout(widths=c(1.3, 1))) /
        (p8b) /
        (p8c + p8d + plot_layout(widths=c(1, 1))) +
        plot_layout(heights=c(1.05, 0.9, 1.25)) +
        plot_annotation(tag_levels="a")

save_fig(fig8, "Figure8_clinical_and_composition", 183, 215)
write.csv(pa_df, file.path(FD, "SourceData_Figure8a.csv"), row.names=FALSE)
write.csv(pb_long, file.path(FD, "SourceData_Figure8b.csv"), row.names=FALSE)
write.csv(pc_long, file.path(FD, "SourceData_Figure8c.csv"), row.names=FALSE)
write.csv(pd_long, file.path(FD, "SourceData_Figure8d.csv"), row.names=FALSE)
write.csv(pe_df, file.path(FD, "SourceData_Figure8e.csv"), row.names=FALSE)

## ---- 8. Build updated Figure 1 overview (Figure1_overview_v6) that includes the 428-biopsy clinical evaluation ----
f1_pan <- readRDS(file.path(FD, "panels_fig1.rds"))
bx <- function(x, y, w, h, lab, fill, col = "black", sz = TX, ff = "plain", tc = "black") list(
  annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill, colour = col, linewidth = .25),
  annotate("text", x = x, y = y, label = lab, size = sz, lineheight = .9, fontface = ff, colour = tc))
hd <- function(x, lab) annotate("text", x = x, y = 10.55, label = lab, size = fs(7), fontface = "bold")
ar <- function(x0, x1, y = 5.2) annotate("segment", x = x0, xend = x1, y = y, yend = y, arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = .45, colour = "grey35")
lt <- function(c) grDevices::adjustcolor(c, alpha.f = .16)
pa_v6 <- ggplot() + coord_cartesian(xlim = c(0, 20), ylim = c(0, 11), expand = FALSE) + theme_void() + theme(plot.tag = element_text(size = 9, face = "bold")) +
  hd(2.6, "Gene-set library (n = 48)") + hd(7.7, "Normal & clinical kidney data") + hd(12.75, "Evaluation & controls") + hd(17.6, "Outputs") +
  bx(2.6, 8.6, 5.0, 2.2, "A  Kidney age-derived (12)\nRodwell, Melk, Rowland, Kainz,\nRAAG, MHC-II, KAS-94/82", lt(cls_col[["A"]]), cls_col[["A"]]) +
  bx(2.6, 5.6, 5.0, 2.2, "B  General ageing (7)\nde Magalhães, Palmer, Sood,\nGenAge, CellAge, SenMayo, Fridman", lt(cls_col[["B"]]), cls_col[["B"]]) +
  bx(2.6, 2.4, 5.0, 2.7, "C  'Ageing-related' disease\nmodels (29): DKD, IgAN, lupus,\ntransplant, RCC; 268 records\nscreened → 29 full gene lists", lt(cls_col[["C"]]), cls_col[["C"]]) +
  bx(7.7, 9.05, 4.7, 1.45, "6 GEO donor cohorts: 441 donors\n(9 series, 498 samples reviewed)", lt(src_col[["G6"]]), src_col[["G6"]]) +
  bx(7.7, 7.25, 4.7, 1.35, "TCGA normal (129); CPTAC\nadjacent (75); GTEx cortex (85)", lt(src_col[["TCGA"]]), src_col[["TCGA"]]) +
  bx(7.7, 5.35, 4.7, 1.65, "Clinical biopsies (n = 428):\nERCB 9 CKD diseases (190),\nGSE175759 eGFR (87),\nGSE147451 graft outcome (151)", lt("#C87028"), "#C87028") +
  bx(7.7, 2.35, 4.7, 2.6, "Leakage rule: a dataset used to\nderive a set is never used to\nevaluate it (KAS: no GEO;\nRowland: no TCGA)", "white", "grey40", ff = "italic") +
  bx(12.75, 8.75, 4.8, 1.9, "Gene-level enrichment vs\n2,000 expression-matched\nrandom sets & camera VIF", "grey96", "grey40") +
  bx(12.75, 6.45, 4.8, 1.9, "Signed score vs age;\nrandom-effects pooling;\n500 random signatures", "grey96", "grey40") +
  bx(12.75, 4.25, 4.8, 1.8, "Gene & score composition\nadjustment + KPMP snRNA-seq\ndeconvolution (MuSiC, Bisque)", "grey96", "grey40") +
  bx(12.75, 2.1, 4.8, 1.8, "Co-primary without CPTAC;\npower by set size;\n47 GTEx tissues", "grey96", "grey40") +
  bx(17.6, 8.3, 4.1, 2.4, "Benchmark scorecard,\npower & tissue source\n(Figs 2–7)", "white", "black", sz = TXL) +
  bx(17.6, 5.3, 4.1, 2.4, "Clinical evaluation in\nCKD & transplantation\n(Fig. 8)", "white", "black", sz = TXL) +
  bx(17.6, 2.3, 4.1, 2.4, "kidneyAgeBench v0.2.0\nR package (Box 1, Fig. S1)", "grey25", "grey25", sz = TXL, ff = "bold", tc = "white") +
  ar(5.1, 5.35) + ar(10.05, 10.35) + ar(15.15, 15.5)
fig1_v6 <- pa_v6 / (f1_pan$pb + f1_pan$pc + plot_layout(widths = c(1, 1.05))) + plot_layout(heights = c(1.05, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig1_v6, "Figure1_overview_v6", 183, 140)
cat("Figure 8 and all revision upgrade tables saved successfully.\n")
