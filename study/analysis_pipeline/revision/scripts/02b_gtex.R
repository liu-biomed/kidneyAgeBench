# Reviewer 1 points 1-2: complete GTEx score-level analysis, arm decomposition, technical covariates,
# single-sample scoring, age-binning sensitivity, and correlation-aware nulls for gene-level concordance
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
set.seed(20261030)
b <- load_bundle(); harm <- b$harm; gt <- b$gtex; P <- readRDS(MP); rob <- P$rob
up <- rob$gene[rob$direction=="up"]; dn <- rob$gene[rob$direction=="down"]
m <- gt$meta; m$vent <- as.integer(m$DTHHRDY==0); m$isch_h <- m$SMTSISCH/60
s <- score_z(gt$expr, up, dn); s$rank <- score_rank(gt$expr, up, dn)
boot_rho <- function(x, y, B=5000) quantile(replicate(B, { i <- sample(length(x), replace=TRUE); cor(x[i], y[i], method="s") }), c(.025,.975))
rows <- list()
for (a in c("score","up","dn","rank")) {
  v <- s[[a]]; sp <- spear(v, m$age); ci <- boot_rho(v, m$age)
  f0 <- lm(scale(v) ~ I(m$age/10) + m$sex); f1 <- lm(scale(v) ~ I(m$age/10) + m$sex + m$vent + m$isch_h + m$SMRIN)
  nv <- m$vent==0; sp_nv <- spear(v[nv], m$age[nv])
  rows[[a]] <- data.frame(axis=a, n=length(v), spearman_rho=sp[1], boot_lo=ci[1], boot_hi=ci[2], p=sp[2],
    beta_per_decade_SD=coef(f0)[2], beta_lo=confint(f0)[2,1], beta_hi=confint(f0)[2,2], p_beta=summary(f0)$coef[2,4],
    beta_adj_tech=coef(f1)[2], p_adj_tech=summary(f1)$coef[2,4],
    n_nonventilator=sum(nv), rho_nonventilator=sp_nv[1], p_nonventilator=sp_nv[2])
}
gscore <- do.call(rbind, rows); print(format(gscore, digits=3), row.names=FALSE); wr(gscore, "R1_GTEx_score_level.csv")
# Which technical factors track the score?
tech <- data.frame(factor=c("ventilator (Hardy 0)","ischaemic time (h)","RIN"),
  rho=c(cor(s$score, m$vent, method="s"), cor(s$score, m$isch_h, method="s"), cor(s$score, m$SMRIN, method="s", use="complete")),
  p=c(spear(s$score, m$vent)[2], spear(s$score, m$isch_h)[2], spear(s$score, m$SMRIN)[2]),
  rho_with_age=c(cor(m$age, m$vent, method="s"), cor(m$age, m$isch_h, method="s"), cor(m$age, m$SMRIN, method="s", use="complete")))
print(tech); wr(tech, "R1_GTEx_technical_factors.csv")

## Age-resolution sensitivity: coarsen discovery ages to GTEx-style decade midpoints and to GTEx age range
mid <- function(a) pmin(pmax(floor(a/10)*10 + 5, 25), 75)
bin <- do.call(rbind, lapply(names(harm), function(g){ h <- harm[[g]]; sc <- score_z(h$expr, up, dn)$score; a <- h$meta$age
  inr <- a >= 20 & a < 80
  data.frame(cohort=g, n=length(a), age_sd=sd(a), rho_exact=cor(sc, a, method="s"), rho_decade_bin=cor(sc, mid(a), method="s"),
             n_20_79=sum(inr), rho_20_79_binned=if (sum(inr) > 5) cor(sc[inr], mid(a[inr]), method="s") else NA) }))
bin <- rbind(bin, data.frame(cohort="GTEx", n=nrow(m), age_sd=sd(m$age), rho_exact=NA, rho_decade_bin=gscore$spearman_rho[1], n_20_79=nrow(m), rho_20_79_binned=NA))
print(format(bin, digits=3), row.names=FALSE); wr(bin, "R1_age_binning_sensitivity.csv")

## Gene-level concordance: binomial vs age-label permutation (preserves gene-gene correlation) vs effective-number binomial
ge <- compute_effects(gt$expr, m, TRUE); rownames(ge) <- ge$gene
g93 <- intersect(rob$gene, ge$gene); disc <- setNames(rob$z, rob$gene)[g93]
obs_conc <- sum(sign(ge[g93,"r"]) == sign(disc)); obs_cor <- cor(ge[g93,"r"], tanh(disc))
B <- 10000; Xg <- gt$expr[g93,]
perm <- t(replicate(B, { a <- m$age; for (sx in unique(m$sex)) { i <- which(m$sex==sx); a[i] <- a[i][sample(length(i))] }
  e <- compute_effects(Xg, m, TRUE, age=a); c(sum(sign(e$r)==sign(disc)), cor(e$r, tanh(disc))) }))
ev <- eigen(cor(t(Xg)), only.values=TRUE)$values; meff <- (sum(ev))^2/sum(ev^2)   # Bretherton effective number
conc <- data.frame(n_genes=length(g93), concordant=obs_conc, pct=100*obs_conc/length(g93),
  p_binomial_naive=binom.test(obs_conc, length(g93), alternative="greater")$p.value,
  M_eff=meff, p_binomial_Meff=pbinom(round(meff*obs_conc/length(g93))-1, round(meff), 0.5, lower.tail=FALSE),
  p_perm_concordance=(1+sum(perm[,1] >= obs_conc))/(B+1), perm_null_mean=mean(perm[,1]), perm_null_95=quantile(perm[,1], .95),
  effect_cor=obs_cor, p_perm_effect_cor=(1+sum(perm[,2] >= obs_cor))/(B+1))
print(t(conc)); wr(conc, "R1_GTEx_concordance_nulls.csv")
# Random-gene null: 93 genes drawn from the meta-analysed universe with matched discovery |r| distribution
res <- P$res; uni <- intersect(res$gene[!res$gene %in% rob$gene], ge$gene); rr <- setNames(res$r, res$gene)
saveRDS(list(score=s, meta=m, ge=ge), file.path(D, "harmonized", paste0("gtex_scores_", KSET, ".rds")))
