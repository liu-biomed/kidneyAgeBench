# Reviewer 1 points 11, 12, 19: within-kidney coupling of the immune-stromal and metabolic programmes,
# and robustness of the canonical-senescence negative result
source(file.path(path.expand("~/Library/Mobile Documents/com~apple~CloudDocs/sandbox/肾内文章/kidney_aging_meta/aging_cell_revision_20261030"), "scripts", "fun.R"))
suppressMessages({library(fgsea); library(msigdbr)})
set.seed(20261030)
b <- load_bundle(); harm <- b$harm; P <- readRDS(MP); rob <- P$rob; res <- P$res
H <- as.data.frame(msigdbr(species="Homo sapiens", collection="H")); hl <- split(H$gene_symbol, H$gs_name)
k94 <- rob$gene
imm_g <- setdiff(unique(unlist(hl[c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_INFLAMMATORY_RESPONSE")])), k94)
met_g <- setdiff(unique(unlist(hl[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_FATTY_ACID_METABOLISM")])), k94)
mz <- function(e, g){ g <- intersect(g, rownames(e)); z <- t(scale(t(e[g,]))); z[!is.finite(z)] <- NA; colMeans(z, na.rm=TRUE) }
dec <- read.csv(file.path(REV, "..", "tables", "deconv_proportions_bisque_class.csv"))
pcor <- function(x, y, C=NULL){ if (is.null(C) || ncol(as.data.frame(C))==0) return(cor(x, y, method="s"))
  C <- as.data.frame(C); cor(residuals(lm(rank(x) ~ ., data=C)), residuals(lm(rank(y) ~ ., data=C))) }
cb <- function(x, y) if (is.null(x)) as.data.frame(y) else cbind(x, y)
## 1. Coupling
rows <- list()
for (g in names(harm)) { h <- harm[[g]]; e <- h$expr; a <- h$meta$age; n <- ncol(e)
  s <- score_z(e, rob$gene[rob$direction=="up"], rob$gene[rob$direction=="down"]); im <- mz(e, imm_g); me <- mz(e, met_g)
  sexc <- if (g %in% NO_SEX) NULL else data.frame(sex=factor(h$meta$sex))
  base <- if (is.null(sexc)) data.frame(age=a) else data.frame(age=a, sex=sexc$sex)
  dd <- dec[match(first_gsm(h$meta$gsm), dec$gsm), c("epithelial","immune","stroma","endothelial")]; hasdec <- all(!is.na(dd))
  comp <- if (hasdec) cbind(base, dd) else NULL
  rows[[g]] <- data.frame(cohort=g, n=n, deconvolved=hasdec,
    arms_raw=cor(s$up, s$dn, method="s"), arms_age=pcor(s$up, s$dn, base), arms_age_comp=if (hasdec) pcor(s$up, s$dn, comp) else NA,
    path_raw=cor(im, me, method="s"), path_age=pcor(im, me, base), path_age_comp=if (hasdec) pcor(im, me, comp) else NA,
    age_imm=pcor(a, im, sexc), age_met=pcor(a, me, sexc), age_met_given_imm=pcor(a, me, cb(sexc, data.frame(imm=im))),
    age_imm_comp=if (hasdec) pcor(a, im, cb(sexc, dd)) else NA, age_met_comp=if (hasdec) pcor(a, me, cb(sexc, dd)) else NA) }
C <- do.call(rbind, rows); wr(C, "R1_arm_coupling_by_cohort.csv")
cols <- setdiff(names(C), c("cohort","n","deconvolved"))
pooledC <- do.call(rbind, lapply(cols, function(v){ ok <- is.finite(C[[v]]); cbind(measure=v, n_samples=sum(C$n[ok]), pool_rho(C[[v]][ok], C$n[ok] - ifelse(grepl("comp", v), 5, 2))) }))
print(format(pooledC, digits=3), row.names=FALSE); wr(pooledC, "R1_arm_coupling_pooled.csv")
## 2. Senescence: expanded gene-set panel, two enrichment frameworks, and direct overlap
M <- as.data.frame(msigdbr(species="Homo sapiens")); sets <- c("SAUL_SEN_MAYO","FRIDMAN_SENESCENCE_UP","REACTOME_CELLULAR_SENESCENCE",
  "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP","WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP","GOBP_CELLULAR_SENESCENCE",
  "KAMMINGA_SENESCENCE","TANG_SENESCENCE_TP53_TARGETS_UP","REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE",
  "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE","GOBP_REPLICATIVE_SENESCENCE","GOBP_STRESS_INDUCED_PREMATURE_SENESCENCE")
sl <- lapply(sets, function(s) unique(M$gene_symbol[M$gs_name==s])); names(sl) <- sets
rk <- res[res$k >= PAR$kmin, ]; st <- sort(setNames(rk$z, rk$gene), decreasing=TRUE)
fg <- as.data.frame(fgseaMultilevel(sl, st, minSize=5, maxSize=1000, eps=0))[, c("pathway","size","NES","pval","padj")]
U <- res$gene; sig136 <- res$gene[res$fdr < 0.05]
sen <- do.call(rbind, lapply(sets, function(s){ G <- intersect(sl[[s]], U)
  # self-contained, sample-level test: per-cohort mean-z set score vs age, pooled across cohorts
  sc <- sapply(names(harm), function(g){ h <- harm[[g]]; v <- mz(h$expr, G); c(rho=cor(v, h$meta$age, method="s"), n=ncol(h$expr)) })
  pr <- pool_rho(sc["rho",], sc["n",])
  ov94 <- intersect(G, k94); ov136 <- intersect(G, sig136)
  data.frame(set=s, measured=length(G), overlap_KAS94=length(ov94), overlap_genes=paste(ov94, collapse=";"),
    p_overlap_KAS94=phyper(length(ov94)-1, length(G), length(U)-length(G), length(k94), lower.tail=FALSE),
    overlap_FDR05=length(ov136), median_member_r=median(res$r[res$gene %in% G]), pct_positive=100*mean(res$z[res$gene %in% G] > 0),
    sample_level_pooled_rho=pr$rho, sample_level_lo=pr$lo, sample_level_hi=pr$hi, sample_level_p=pr$p, sample_level_I2=pr$I2) }))
sen <- merge(sen, setNames(fg, c("set","fgsea_size","fgsea_NES","fgsea_p","fgsea_padj")), by="set", all.x=TRUE)
sen$sample_level_fdr <- p.adjust(sen$sample_level_p, "BH")
print(format(sen[, -4], digits=2), row.names=FALSE); wr(sen, "R1_senescence_expanded.csv")
# Reference: the same sample-level test for the KAS-94 immune arm and for the immune Hallmark score
ref <- sapply(names(harm), function(g){ h <- harm[[g]]; c(cor(mz(h$expr, imm_g), h$meta$age, method="s"), ncol(h$expr)) })
print(pool_rho(ref[1,], ref[2,]))
