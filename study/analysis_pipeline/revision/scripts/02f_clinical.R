# Reviewer 1 points 20-23: clinical specificity against matched null gene sets and comparator signatures,
# arm-level contributions, and adjustment for generic inflammation/injury programmes
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
suppressMessages({library(msigdbr); library(pROC)})
set.seed(20261030)
cl <- readRDS(file.path(D,"harmonized","clinical_bundle.rds")); P <- readRDS(MP); rob <- P$rob
up <- rob$gene[rob$direction=="up"]; dn <- rob$gene[rob$direction=="down"]
H <- as.data.frame(msigdbr(species="Homo sapiens", collection="H")); hl <- split(H$gene_symbol, H$gs_name)
# Comparators: generic inflammation, IFN-gamma, EMT/fibrosis, OXPHOS (inverted), a tubular-injury panel
injury <- c("HAVCR1","LCN2","SPP1","VCAM1","CCL2","TIMP1","CD44","SOX9","KRT8","KRT18","CLU","VIM","ICAM1","C3","CXCL10")
comp <- list(inflammatory=list(up=hl$HALLMARK_INFLAMMATORY_RESPONSE, dn=character()), ifng=list(up=hl$HALLMARK_INTERFERON_GAMMA_RESPONSE, dn=character()),
             emt=list(up=hl$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, dn=character()), oxphos_inv=list(up=character(), dn=hl$HALLMARK_OXIDATIVE_PHOSPHORYLATION),
             injury=list(up=injury, dn=character()), KAS94=list(up=up, dn=dn), KAS94_up=list(up=up, dn=character()), KAS94_down=list(up=character(), dn=dn))
sc1 <- function(e, s){ u <- intersect(s$up, rownames(e)); d <- intersect(s$dn, rownames(e)); z <- t(scale(t(e[c(u,d),,drop=FALSE]))); z[!is.finite(z)] <- 0
  (if (length(u)) colMeans(z[u,,drop=FALSE]) else 0) - (if (length(d)) colMeans(z[d,,drop=FALSE]) else 0) }
auc1 <- function(y, x) as.numeric(auc(roc(y, x, direction="<", quiet=TRUE)))
## 1. Comparator signatures
cmp <- do.call(rbind, lapply(names(cl), function(co){ e <- cl[[co]]$expr; m <- cl[[co]]$meta; y <- m$grp=="Disease"
  do.call(rbind, lapply(names(comp), function(nm){ s <- sc1(e, comp[[nm]]); r <- roc(y, s, direction="<", quiet=TRUE); ci <- ci.auc(r)
    data.frame(cohort=co, signature=nm, genes=length(intersect(c(comp[[nm]]$up, comp[[nm]]$dn), rownames(e))), AUC=as.numeric(auc(r)), AUC_lo=ci[1], AUC_hi=ci[3],
      egfr_rho=if (co=="GSE175759") cor(s, m$egfr, method="s") else NA, egfr_p=if (co=="GSE175759") spear(s, m$egfr)[2] else NA) })) }))
print(format(cmp, digits=3), row.names=FALSE); wr(cmp, "R1_clinical_comparator_signatures.csv")
## 2. Null gene sets matched on mean expression and variance (per cohort), 1,000 draws
mnull <- do.call(rbind, lapply(names(cl), function(co){ e <- cl[[co]]$expr; m <- cl[[co]]$meta; y <- m$grp=="Disease"
  mu <- rowMeans(e); sdv <- apply(e, 1, sd); bins <- interaction(cut(rank(mu), 10), cut(rank(sdv), 10)); names(bins) <- rownames(e)
  pick <- function(g){ g <- intersect(g, rownames(e)); unlist(lapply(g, function(x){ pool <- setdiff(names(bins)[bins==bins[x]], c(up,dn)); sample(pool, 1) })) }
  obs <- sc1(e, comp$KAS94); oa <- auc1(y, obs); oe <- if (co=="GSE175759") cor(obs, m$egfr, method="s") else NA
  nl <- t(replicate(1000, { s <- sc1(e, list(up=pick(up), dn=pick(dn))); c(auc1(y, s), if (co=="GSE175759") cor(s, m$egfr, method="s") else NA) }))
  data.frame(cohort=co, observed_AUC=oa, null_AUC_median=median(nl[,1]), null_AUC_95=quantile(nl[,1], .95), p_AUC=(1+sum(nl[,1] >= oa))/1001,
             observed_egfr_rho=oe, null_rho_median=median(nl[,2]), null_rho_05=quantile(nl[,2], .05, na.rm=TRUE),
             p_egfr=if (co=="GSE175759") (1+sum(nl[,2] <= oe))/1001 else NA) }))
print(format(mnull, digits=3), row.names=FALSE); wr(mnull, "R1_clinical_matched_null.csv")
## 3. Incremental information: combined vs arms, and KAS-94 after adjustment for generic inflammation/injury scores
inc <- do.call(rbind, lapply(names(cl), function(co){ e <- cl[[co]]$expr; m <- cl[[co]]$meta; y <- as.integer(m$grp=="Disease")
  d <- data.frame(y=y, up=scale(sc1(e, comp$KAS94_up))[,1], dn=scale(sc1(e, comp$KAS94_down))[,1], k=scale(sc1(e, comp$KAS94))[,1],
                  infl=scale(sc1(e, comp$inflammatory))[,1], inj=scale(sc1(e, comp$injury))[,1])
  lr <- function(f0, f1) anova(glm(f0, d, family=binomial), glm(f1, d, family=binomial), test="LRT")[2, "Pr(>Chi)"]
  out <- data.frame(cohort=co, p_down_beyond_up=lr(y~up, y~up+dn), p_up_beyond_down=lr(y~dn, y~up+dn),
             p_KAS94_beyond_inflammation_injury=lr(y~infl+inj, y~infl+inj+k), cor_KAS94_inflammation=cor(d$k, d$infl, method="s"), cor_KAS94_injury=cor(d$k, d$inj, method="s"))
  if (co=="GSE175759") { eg <- m$egfr; out$egfr_partial_rho_adj_infl_inj <- cor(residuals(lm(rank(d$k) ~ d$infl + d$inj)), residuals(lm(rank(eg) ~ d$infl + d$inj)))
    f <- lm(scale(eg) ~ d$k + d$infl + d$inj); out$egfr_beta_KAS94_adj <- coef(f)[2]; out$egfr_p_KAS94_adj <- summary(f)$coef[2,4] }
  else { out$egfr_partial_rho_adj_infl_inj <- NA; out$egfr_beta_KAS94_adj <- NA; out$egfr_p_KAS94_adj <- NA }
  out }))
print(format(inc, digits=3), row.names=FALSE); wr(inc, "R1_clinical_incremental.csv")
## 4. GSE66494 reference material
r66 <- data.frame(cl$GSE66494$meta[, c("gsm","grp","source")], cl$GSE66494$sample_chars[, 2]); names(r66)[4] <- "sample_type"
wr(as.data.frame(table(r66$grp, r66$sample_type)), "R1_GSE66494_reference_material.csv")
