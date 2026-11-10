# Addendum 2, D1-D3 (CPTAC discovery, paired NAT RNA + protein) and UMich technical replication
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "10_cptac_helpers.R"))
suppressMessages({library(fgsea); library(msigdbr)}); set.seed(20261030)
V <- path.expand("~/kidney_aging_data/validation/cptac"); mp <- cptac_rd(V, "bcm-ccrcc-mapping"); e2s <- setNames(mp$gene_name, sub("\\..*", "", mp$gene)); C <- cptac_clin(V)
P <- cptac_sym(cptac_rd(V, "bcm-ccrcc-proteomics.*Normal"), e2s); P <- P[rowMeans(is.finite(P)) >= 0.5, ]
tr <- cptac_rd(V, "bcm-ccrcc-transcriptomics"); tr <- cbind(gene=tr[[1]], tr[, grep("_A$", names(tr))]); R <- cptac_sym(tr, e2s); colnames(R) <- sub("_A$", "", colnames(R)); R <- R[rowMeans(R > 0) >= 0.2, ]
cs <- Reduce(intersect, list(colnames(P), colnames(R), C$case[is.finite(C$age) & C$sex %in% c("Male","Female")])); P <- P[, cs]; R <- R[, cs]; m <- C[match(cs, C$case), ]
cat("paired NAT patients:", length(cs), "| proteins:", nrow(P), "| RNA genes:", nrow(R), "\n")
FC <- full_cov(P, m$sex); sm <- LOCK$outcome$SenMayo
## D1 score level
d1 <- rbind(cbind(test="protein SenMayo ~ age | full model", pspear(m$age, mz(P, sm), FC)),
            cbind(test="protein SenMayo ~ age | full model + RNA SenMayo", pspear(m$age, mz(P, sm), cbind(FC, rnaS=mz(R, sm)))),
            cbind(test="RNA SenMayo ~ age | full model (protein proxies)", pspear(m$age, mz(R, sm), FC)))
print(format(d1, digits=3), row.names=FALSE)
## Gene-level PTAE for all genes with both layers (used for D1 gene level and D2)
gb <- intersect(rownames(P), rownames(R)); cat("genes with both layers:", length(gb), "\n")
Zc <- residuals(lm(rank(m$age) ~ ., data=FC))
ptae <- do.call(rbind, lapply(gb, function(g){ p <- P[g,]; r <- R[g,]; ok <- is.finite(p) & is.finite(r) & complete.cases(FC); if (sum(ok) < 30) return(NULL)
  Cg <- cbind(FC[ok,,drop=FALSE], mrna=rank(r[ok])); a <- residuals(lm(rank(m$age[ok]) ~ ., Cg)); pp <- residuals(lm(rank(p[ok]) ~ ., Cg))
  a0 <- residuals(lm(rank(m$age[ok]) ~ ., FC[ok,,drop=FALSE])); p0 <- residuals(lm(rank(p[ok]) ~ ., FC[ok,,drop=FALSE])); r0 <- residuals(lm(rank(r[ok]) ~ ., FC[ok,,drop=FALSE]))
  data.frame(gene=g, n=sum(ok), ptae=cor(a, pp), prot_age=cor(a0, p0), rna_age=cor(a0, r0), rna_prot=cor(rank(p[ok]), rank(r[ok]))) }))
write.csv(ptae, file.path(O, "CPTAC_PTAE_by_gene.csv"), row.names=FALSE)
sg <- ptae[ptae$gene %in% sm, ]
d1g <- data.frame(members=nrow(sg), positive_ptae=sum(sg$ptae > 0), sign_test_p=binom.test(sum(sg$ptae > 0), nrow(sg))$p.value,
                  median_ptae=median(sg$ptae), median_prot_age=median(sg$prot_age), median_rna_age=median(sg$rna_age),
                  genome_median_ptae=median(ptae$ptae), wilcox_vs_genome_p=wilcox.test(sg$ptae, ptae$ptae[!ptae$gene %in% sm])$p.value)
print(d1g); write.csv(d1, file.path(O, "D1_score_level.csv"), row.names=FALSE); write.csv(d1g, file.path(O, "D1_gene_level.csv"), row.names=FALSE)
## D2 enrichment of PTAE
c2 <- as.data.frame(msigdbr(species="Homo sapiens", collection="C2")); hh <- as.data.frame(msigdbr(species="Homo sapiens", collection="H"))
prim <- list(NABA_SECRETED_FACTORS=unique(c2$gene_symbol[c2$gs_name=="NABA_SECRETED_FACTORS"]), NABA_CORE_MATRISOME=unique(c2$gene_symbol[c2$gs_name=="NABA_CORE_MATRISOME"]))
sec <- list(SAUL_SEN_MAYO=sm, HALLMARK_OXIDATIVE_PHOSPHORYLATION=LOCK$outcome$OXPHOS)
st <- sort(setNames(ptae$ptae, ptae$gene), decreasing=TRUE)
f1 <- as.data.frame(fgsea(prim, st, minSize=10, maxSize=2000, nPermSimple=10000)); f1$tier <- "primary"; f1$padj <- p.adjust(f1$pval, "BH")
f2 <- as.data.frame(fgsea(sec, st, minSize=10, maxSize=2000, nPermSimple=10000)); f2$tier <- "secondary"
f3 <- as.data.frame(fgsea(split(hh$gene_symbol, hh$gs_name), st, minSize=10, maxSize=2000, nPermSimple=10000)); f3$tier <- "exploratory_hallmark"
F <- rbind(f1, f2, f3)[, c("tier","pathway","size","NES","pval","padj")]; F <- F[order(F$tier, F$pval),]
print(format(head(F[F$tier!="exploratory_hallmark",], 10), digits=3), row.names=FALSE); cat("top Hallmark (exploratory):\n"); print(format(head(F[F$tier=="exploratory_hallmark",], 12), digits=3), row.names=FALSE)
write.csv(F, file.path(O, "D2_PTAE_enrichment.csv"), row.names=FALSE)
## D3 decoupling with age (median split; age-permutation null)
young <- m$age <= median(m$age)
coup <- function(yg){ zy <- sapply(gb, function(g) suppressWarnings(cor(P[g, yg], R[g, yg], method="s", use="complete")))
  zo <- sapply(gb, function(g) suppressWarnings(cor(P[g, !yg], R[g, !yg], method="s", use="complete"))); d <- atanh(pmin(pmax(zo, -.999), .999)) - atanh(pmin(pmax(zy, -.999), .999))
  c(median_diff=median(d, na.rm=TRUE), young_median=median(zy, na.rm=TRUE), old_median=median(zo, na.rm=TRUE)) }
obs <- coup(young); nul <- replicate(2000, coup(sample(young))["median_diff"])
d3 <- data.frame(n_young=sum(young), n_old=sum(!young), age_split=median(m$age), young_median_r=obs["young_median"], old_median_r=obs["old_median"],
                 median_fisher_diff_old_minus_young=obs["median_diff"], perm_p_two_sided=(1+sum(abs(nul) >= abs(obs["median_diff"])))/2001)
print(d3); write.csv(d3, file.path(O, "D3_decoupling.csv"), row.names=FALSE)
